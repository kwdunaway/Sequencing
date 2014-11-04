#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 11/3/2014
#
# This script checks to see what is in each tab
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Input file name
" unless @ARGV == 1;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";


###################
# Another section #
###################

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    print "NEW LINE\n";
    for(my $n = 0; $n < @line; $n++){
    	print "tab $n" , ":\t" , $line[$n] , "\n";
    }
}
close IN;
