#!/usr/bin/perl

my $fastaFile = shift;

my @Seq = ();
my $total = 0;
my @Name = ();
my @Abund = ();

open(FILE, $fastaFile) or die "Can't open $fastaFile\n";

my $seq = "";

while($line = <FILE>){ 
    chomp($line);
    
    if($line =~ />(.*)/){
	
	@tokens = split(/_/,$line);
	$size = scalar(@tokens);

	$Abund[$count] = $tokens[$size-1];
	$Abund[$count] =~ s/[^0-9]//g;
	pop(@tokens);
	$Name[$count] = join("_", @tokens);
	
	if($seq ne ""){
	    $Seq[$count - 1] = $seq;

	    $seq = "";
	}

	$count++;
    }
    else{
	$seq .= $line;
    }
}

$Seq[$count - 1] = $seq;

$total = $count;

for($i = 0; $i < $total; $i++){
	if($Abund[$i] > 0){
		printf STDOUT "$Name[$i];size=$Abund[$i];\n$Seq[$i]\n";
	}
}
close(FILE);
