#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: filterIDIC.pl
# Version: 1.0
# Last Updated: 3/31/2013
#
# This script filters an IDIC file for any -- on single lines that screws up stuff.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Infile
    2) Outfile
" unless @ARGV == 2;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


#############
# Main Loop #
#############

while (<IN>)
{
    chomp;
    if($_ eq "--"){}
    else { print OUT $_  , "\n";} 
}

