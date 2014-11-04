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
    1) Input file name
    2) Output file name
" unless @ARGV == 2;

my $fastqinfile = shift(@ARGV);
open(FASTQ, "<$fastqinfile") or die "cannot open $fastqinfile fastqinfile";
my $searchinfile = shift(@ARGV);
open(SEARCH, "<$searchinfile") or die "cannot open $searchinfile searchinfile";

######################
# Load Search String #
######################

my %SearchHash;
my @search;
while(<SEARCH>){
    chomp;
    $SearchHash{$_} = 0;
    push(@search;$_);
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
    for($n = 0; $n < @search; $n++){
    	if($seq =~ /$search[$n]/){$SearchHash{$search[$n]}++;}
    }
}

########################
# Print Search results #
########################

for($n = 0; $n < @search; $n++){
    print $search[$n] , "\t" , $SearchHash{$search[$n]} , "\n";
}
