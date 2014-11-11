#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 11/4/2014
#
# This script takes in a file and extracts a subset of a column in the file
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Output prefix (ends with found string.txt)
    2) Input FASTQ file
    3) Input SEARCH file
" unless @ARGV == 3;

my $outprefix = shift(@ARGV);
my $fastqinfile = shift(@ARGV);
open(FASTQ, "<$fastqinfile") or die "cannot open $fastqinfile fastqinfile";
my $searchinfile = shift(@ARGV);
open(SEARCH, "<$searchinfile") or die "cannot open $searchinfile searchinfile";

######################
# Load Search String #
######################

my %SearchHash;
my %StringHash;
my @search;
while(<SEARCH>){
    chomp;
    my $rev = reverse($_);
    my $revcomp = complement($rev);
    my $comp = complement($_);

    push(@search,$_);
    push(@search,$rev);
    push(@search,$revcomp);
    push(@search,$comp);
    $SearchHash{$_} = 0;
    $SearchHash{$rev} = 0;
    $SearchHash{$revcomp} = 0;
    $SearchHash{$comp} = 0;
    
}

#####################
# Search FASTQ file #
#####################

while (<FASTQ>)
{
    my $seq = <FASTQ>;
    chop($seq);
    <FASTQ>;
    <FASTQ>;
    for(my $n = 0; $n < @search; $n++){
    	if($seq =~ /$search[$n]/){
    		$SearchHash{$search[$n]}++;
    		if($SearchHash{$search[$n]} % 10000 == 0){
	    		print $search[$n] , "\t" , $SearchHash{$search[$n]}, "\n";
	    		$StringHash{$search[$n]}{$SearchHash{$search[$n]}} = $seq;
	    	}
    	}
    }
}

########################
# Print Search results #
########################

for(my $n = 0; $n < @search; $n++){
    print $search[$n] , "\t" , $SearchHash{$search[$n]} , "\n";
#    if($SearchHash{$search[$n]} > 1000){
#    	my $outfile = $outprefix . "_" . $search[$n] . ".txt";
#		open(OUT, ">$outfile") or die "cannot open $outfile outfile";    	
#		for(my $i = 0; $i < $SearchHash{$search[$n]}; $i++){
#			print OUT $StringHash{$search[$n]}{$SearchHash{$search[$n]}} , "\n";
#		}
#	}
}


###############
# Subroutines #
###############

sub complement {
    my $dna = shift;
    $dna =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $dna;
}

