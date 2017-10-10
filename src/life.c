//
//
//  Mark Sutherland - 997332989
//  David Biancolin - 997516512
//
//

/****************************************************************************************
 * ALGORITHM FOR OPTIMIZED GAME OF LIFE ROUTINE RESIDES HERE.
 *         
 * Optimization list:
 *  - Changed array access pattern from row-major to column-major order because of 
 *    the way that the board is set in memory, and determined by BOARD macro.
 *
 *  - Used pthreads library (4 threads total) to perform the multiple generation 
 *    calculations faster. Synchronized generational time using shared barrier.
 *
 *  - Entirely re-worked the process of looking up a particular iteration's neighbor
 *    count and computing its alive/dead status. Uses shared lookup table, indexed by
 *    the particular combination of neighbors, to determine next generational alive
 *    or dead status. INDEXED AS FOLLOWS:
 *          Suppose a 3x3 section of the board is as follows: A = alive, D = dead
 *
 *          D ---- D ---- A
 *          D ---- D ---- D
 *          A ---- A ---- D
 *
 *          This will be represented in the lookup table by the following 1 byte 
 *          value: 0 0 1 0 0 0 1 1 0. As such, looking up this particular key allows
 *          us to quickly determine in constant time what the status of this particular
 *          cell will be (in this case alive since it has exactly 3 neighbors).
 *
 *  - Manually unrolled the computation of new cells by 4x (by realizing that each 
 *    side-by-side pair of cells shares 3 pairs of computational values needed). Taking
 *    the same example above, the pair of elements <top-middle,top_right> is needed for 
 *    the computation of both the current centre cell and what is shown as the middle-right
 *    cell. Saved reads by taking these value pairs and bit-shifting them into variables "key1" 
 *    and "key2-4". 
 *
 *  - When moving down columns, each new row only requires 2 more memory accesses/loads
 *    per element (amortized over the 4x unroll, 10 reads raw but 1/2 the number of columns
 *    to compute). This is done by shifting the 3 MSB's of the lookup key out, and then using
 *    more bit-shifts to bring in the new values that are needed in the 3 LSB's.
 *
 *  - Removed costly "mod" check for the majority of cell computation. Since this is only 
 *    required when we are at the boundary of the board, we can do 1 check for this condition
 *    at the beginning of each column/row pair and save the repeated calls to mod.
 *******************************************************************************************/

#include "life.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>
#define SWAP_BOARDS( b1, b2 )  do { \
  char* temp = b1; \
  b1 = b2; \
  b2 = temp; \
} while(0)

#define BOARD( __board, __i, __j )  (__board[(__i) + LDA*(__j)])

#define KEY_LENGTH 9 

struct tdata {
    int cindex;
    int colsToUse;
    int nrows;
    int generations;
    pthread_barrier_t* generation_barrier;
    char* lookup;
    char* inboard;
    char* outboard;
};

/*****************************************************************************
 * Helper function definitions
 ****************************************************************************/
void *thread_life(void *args)
{

    /* struct pointer is passed in as arg - contains relevant data that all of 
     * the threads need (such as sizes, board pointers, the shared lookup table,
     * etc etc).
     */
    struct tdata* my_data = (struct tdata*)args;
    char* lookup = my_data -> lookup;
    char* inboard = my_data -> inboard;
    char* outboard = my_data -> outboard;
    const int nrows = my_data->nrows;
    const int ncols = my_data->nrows;
    const int generations = my_data->generations;

    pthread_barrier_t * generation_barrier = my_data->generation_barrier;

    const int LDA = nrows;
    int i,j,k;
    int endCol = my_data->colsToUse;

    for(k=0;k<generations;k++) { // for each generation do the same quartile of the board
        for(j = my_data->cindex ; j<endCol;j+=4) {

            //Essentially we are striping a column of the input matrix. With a parallel unroll
            //of 4. 4 keys are generated using the scheme above to look up the result for a single
            //element. 
            i=0;
            int key_1 = 0;
            int key_2 = 0;
            int key_3 = 0;
            int key_4 = 0;
            //These define the bounds of our stripe. 
            int jwest = j-1;
            int jeast = j+1; 
            int jeastmost=j+4;

            //Replace calls to mod with these checks; called once per stripe. 
            if (j== 0)
                jwest = ncols-1; 

            if (j >= ncols-4) 
                jeastmost = 0;
                

            int inorth = nrows-1; 
            //Load the first row of values into the keys (3 bits each, bits 8 - 6)
            int shared = (BOARD (inboard, inorth, j)<<7) + (BOARD (inboard, inorth, jeast)<<6); 
            key_1 = (BOARD (inboard, inorth, jwest)<<8) + shared;
            shared = (shared << 1) + (BOARD (inboard, inorth, jeast+1)<<6);
            key_2 = shared;
            shared = (shared << 1) + (BOARD(inboard, inorth, jeast+2)<<6); 
            key_3 = shared & 511; 
            shared = (shared << 1) + (BOARD(inboard, inorth, jeastmost)<<6);
            key_4 = shared & 511; 

            //Load the second row of values bits 5 through 3
            shared = (BOARD (inboard, i, j)<<4) + (BOARD (inboard, i, jeast)<<3); 
            key_1 += (BOARD (inboard, i, jwest)<<5) + (shared);
            shared = (shared << 1) + (BOARD (inboard, 0, jeast+1)<<3);
            key_2 += shared & 63;
            shared = (shared << 1) + (BOARD(inboard, 0, jeast+2)<<3); 
            key_3 += shared & 63; 
            shared = (shared << 1) + (BOARD(inboard, 0, jeastmost)<<3);
            key_4 += shared & 63; 

            //Load third row of values (bits 2 through 0) 
            shared = (BOARD (inboard, 1, j)<<1) + BOARD (inboard, 1, jeast); 
            key_1 += (BOARD (inboard, 1, jwest)<<2) + shared; 
            shared = (shared << 1) + BOARD(inboard, 1, jeast+1);
            key_2 += shared;
            shared = (shared << 1) + BOARD(inboard, 1, jeast+2); 
            key_3 += shared & 7; 
            shared = (shared << 1) + BOARD(inboard, 1, jeastmost);
            key_4 += shared & 7;
            //Perform the alivep look ups with the generated keys
            BOARD(outboard, 0, j) = lookup[key_1];
            BOARD(outboard, 0, jeast) = lookup[key_2];
            BOARD(outboard, 0, jeast+1) = lookup[key_3];
            BOARD(outboard, 0, jeast+2) = lookup[key_4];

            //Note this loop omits the final row of the array
            for (i = 1; i < nrows-1; i+=1) 
            {
                int isouth =i+1;
                //Shift one row down. Push the northern most alive bits out of the keys with 
                //a bitshift of 3 and mask. 
                key_1 = (key_1 << 3); 
                key_2 = (key_2 << 3); 
                key_3 = (key_3 << 3); 
                key_4 = (key_4 << 3);

                //Do as above, load the southern most values (newest)
                int shared = ((BOARD (inboard, isouth, j))<<1) + BOARD (inboard, isouth, jeast);
                key_1 += (BOARD (inboard, isouth, jwest)<<2) + shared; 
                shared = (shared << 1) + BOARD(inboard, isouth, jeast+1);
                key_2 += shared;
                shared = (shared << 1) + BOARD(inboard, isouth, jeast+2); 
                key_3 += (shared & 7); 
                shared = (shared << 1) + BOARD(inboard, isouth, jeastmost);
                key_4 += (shared & 7); 

                //Lookup and set values. 
                BOARD(outboard, i, j) = lookup[key_1&511];
                BOARD(outboard, i, jeast) = lookup[key_2&511];
                BOARD(outboard, i, jeast+1) = lookup[key_3&511];
                BOARD(outboard, i, jeast+2) = lookup[key_4&511];

            }

            //Handle the final row where south wraps around to 0 
            key_1 = ((key_1 << 3)&511); 
            key_2 = ((key_2 << 3)&511);
            key_3 = ((key_3 << 3)&511);
            key_4 = ((key_4 << 3)&511);

            shared = ((BOARD (inboard, 0, j))<<1) + BOARD (inboard, 0, jeast);
            key_1 += (BOARD (inboard, 0, jwest)<<2) + shared; 
            shared = (shared << 1) + BOARD(inboard, 0, jeast+1);
            key_2 += shared;
            shared = (shared << 1) + BOARD(inboard, 0, jeast+2); 
            key_3 += (shared & 7); 
            shared = (shared << 1) + BOARD(inboard, 0, jeastmost);
            key_4 += (shared & 7); 
        
            BOARD(outboard, i, j) = lookup[key_1];
            BOARD(outboard, i, jeast) = lookup[key_2];
            BOARD(outboard, i, jeast+1) = lookup[key_3];
            BOARD(outboard, i, jeast+2) = lookup[key_4];
        }        
        SWAP_BOARDS(inboard,outboard); // this works because these are local scope variables

        /* Now wait for all threads to finish this generation. Will be re-woken when the barrier hits 4. */
        pthread_barrier_wait(generation_barrier);
    }
    pthread_exit(NULL);
}


/*****************************************************************************
 * Final parellelized game of life implementation. 
 ****************************************************************************/
char*
game_of_life (char* outboard, 
	      char* inboard,
	      const int nrows,
	      const int ncols,
	      const int gens_max)
{
    int i, j,rc;

    /* if world size < 32, don't run in parallel. */
    if(nrows < 32)
        return sequential_game_of_life(outboard,inboard,nrows,ncols,gens_max);

    /* else proceed with massively parallel implementation */

    /* Create and initialize barrier. */
    pthread_barrier_t * generation_barrier = malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(generation_barrier,NULL,4);

    char *lookup = malloc(512*sizeof(char)); 
    struct tdata *TDATA = malloc(4*sizeof(struct tdata)); // used to pass various args to the forked threads.
    pthread_t threads[4]; 

    /* Creates the lookup table for the various indices corresponding to 3x3 arrays on the board. */
    for (i=0; i<512; i++) { 
        char alive = 0; 
        char neighbor_count = 0; 
        for (j=0; j<KEY_LENGTH; j++) { // for all possible shifts of this index value
            if ((i>>j)&1){
                if (j == 4) 
                    alive = 1; 
                else
                    neighbor_count++; 
             }
        }

        lookup[i] = alivep(neighbor_count,alive);
     }

    /* Sets up the structs that we are going to pass to the threads. */

    int quartile = nrows/4; /* Yes we are doing this because we can assume the input is a power of 2 */
    for(i=0;i<4;i++) {
        TDATA[i].inboard = inboard;
        TDATA[i].outboard = outboard;
        TDATA[i].cindex = i * quartile;
        TDATA[i].colsToUse = (i+1) * quartile;
        TDATA[i].generations = gens_max;
        TDATA[i].generation_barrier = generation_barrier;
        TDATA[i].lookup = lookup;
        TDATA[i].nrows = nrows;
    }

    /* Create all of the worker threads, and wait on all of them to exit 
     * which will happen after every generation is complete. */
    for(i=0;i<4;i++) {
        rc = pthread_create(&threads[i],NULL,thread_life,(void*)&TDATA[i]);
        if(rc) {
            printf("ERR: pthread_create() for worker thread # %d returned: %d\n",i,rc);
            exit(-1);
        }
    }
    
    /* Wait on all the workers to finish. They are auto-synchronized by using the
     * thread barrier. */
    for(i=0;i<4;i++) 
        pthread_join(threads[i],NULL);

    /* Cleanup. */
    pthread_barrier_destroy(generation_barrier);
    free(lookup);
    free(TDATA);

    /* Make sure to check if we should return the inboard or the outboard pointer!!!
     * Since the threads are swapping their local pointers, if the number of gens is an odd number
     * return the outboard ptr. Even number means return the inboard. */
    if (gens_max % 2) // odd
        return outboard;
    else
        return inboard;
}
