#!/usr/bin/perl

my $amp = shift;
my $full = shift;
my $marginL = shift;
my $marginR = shift;
my @name;
my @patterns;
my $count = 0;

open(FILE, $amp) or die "Can't open $amp\n";

while($line = <FILE>){ 
	chomp($line);
	if($line =~ />(.*)/){
		$name[$count] = $line;
	}
	else{
		$patterns[$count] = $line;
		$count++;
	}
}

close(FILE);

open(FILE, $full) or die "Can't open $full\n";

$count = 0;

while($line = <FILE>){ 
	chomp($line);
	@tokens = split(" ",$line);
	if($name[$count] eq $tokens[0]){
		printf STDOUT "$name[$count]\n";
	}
	elsif($line =~ /(.{$marginL})($patterns[$count])(.{$marginR})/g ){
      	my $left_margin     = $1;
    	my $result          = $2;
       	my $right_margin    = $3;
        		
       	printf STDOUT "$left_margin$result$right_margin\n";
       	$count++;
   	}
}

close(FILE);