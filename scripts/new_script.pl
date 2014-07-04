#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: 
# Version: 0.1
# Last Updated: 
#
# <brief description of script's function>
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) 
    2) 
    3) 
" unless @ARGV == 3;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";


###################
# Another section #
###################

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
}

# 2D arrays
@{$OutputArray[$linenum]} = @line;
