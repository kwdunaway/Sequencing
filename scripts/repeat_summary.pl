#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 9-25-2014
#
# This script will take the output of an overlap with repeats and summarize the results.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Output file
    2+) Input file(s)
" unless @ARGV > 1;

my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
print OUT "Class";
my @Infiles = @ARGV;

#Global Variables
my %Class;
my %Family;

###################
# Another section #
###################

for(my $n = 0; $n < @Infiles; $n++){
	my $infile = $Infiles[$n];
	open(IN, "<$infile") or die "cannot open $infile infile";
	while (<IN>)
	{
	    chomp;
	    my @line = split ("\t", $_);
	    if(defined $Class{$line[1]}{$infile}) {$Class{$line[1]}{$infile}++;}
	    else{$Class{$line[1]}{$infile} = 1;}
	    if(defined $Family{$line[2]}{$infile}) {$Family{$line[2]}{$infile}++;}
	    else{$Family{$line[2]}{$infile} = 1;}
	}
	print OUT "\t" , $infile;
}


#print OUT "Class" , "\t" , "Count" , "\n";
print OUT "\n";
foreach my $key (sort keys %Class){
	print OUT $key;
	for(my $n = 0; $n < @Infiles; $n++){
		my $infile = $Infiles[$n];
	    if(defined $Class{$key}{$infile}) {print OUT "\t" , $Class{$key}{$infile};}
	    else {print OUT "\t" , 0;}
	}
	print OUT "\n";
}

print OUT "\nFamily\n";
foreach my $key (sort keys %Family){
	print OUT $key;
	for(my $n = 0; $n < @Infiles; $n++){
		my $infile = $Infiles[$n];
	    if(defined $Family{$key}{$infile}) {print OUT "\t" , $Family{$key}{$infile};}
	    else {print OUT "\t" , 0;}
	}
	print OUT "\n";
}
