open (FH, "time.out") 
	or die "WTF";

my $average = 0; 

while (<FH>) {
	if ($_ =~ m/([0-9.]*)elapsed/) {
		$average = $average + $1; 
		}
	} 

$average = $average/5;

print "$average\n";
close $fh


