#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: Excel_to_bed.pl
# Version: 1.0
# Last Updated: 
#
# Takes a tab separated excel sheet (saved as a .txt file), pulls out the chromosomal
# positions, and creates a bed file from them.  The chr positions need to be in
# chr##:(start)-(end) format in order for this to work.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input file name
    2) Output file name
    3) Column number of chr position
    4) Track name
" unless @ARGV == 4;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";
my $trackname = shift(@ARGV);
my $color = shift(@ARGV);

my $header = "track name=\"" . $trackname . "\" description=\"" . $trackname . "\" itemRgb=\"On\"";
print OUT $header;

###################
# Another section #
###################

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    my @chr = split (":", $line[1]);
    my @pos = split("-", $chr[1]);

	# Print format
	# chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
    print OUT "\n", $chr[0] , "\t" , $pos[0], "\t" , $pos[1], "\t", $line[0], "\t0\t+\t", $pos[0], "\t" , $pos[1], "\t", $color;
}


