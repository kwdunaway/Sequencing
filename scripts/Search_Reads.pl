#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: Search_Reads.pl
# Version: 0.1
# Last Updated: 2/21/2014
#
# This script takes a fastq file and searches the reads for a particular sequence, then
# prints any matches to the screen.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Fastq file (.fq)
    2) String to match
    3) Line to match it to (1-4)
" unless @ARGV == 3;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $matchstring = shift(@ARGV);
my $linenumber = shift(@ARGV);

#my $output_filename = shift(@ARGV);
#open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


#############
# Main Loop #
#############

while (<IN>)
{
    my @line;
    $line[1] = $_;
    $line[2] = <IN>;
    $line[3] = <IN>;
    $line[4] = <IN>;
    
    if ($line[$linenumber] =~ /$matchstring/){
    	print $line[1] . $line[2] . $line[3] . $line[4];
    }
}

