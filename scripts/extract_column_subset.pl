#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 11/3/2014
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
    3) column
    4) substring number (ex: 5 is first 5 chars, ex2: -10 is last 10 chars)
" unless @ARGV == 4;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
my $column = shift(@ARGV);
my $lastnum = shift(@ARGV);


###################
# Another section #
###################

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    my $stringout = substr($line[$column],$lastnum);
    print OUT $stringout , "\n";
}

