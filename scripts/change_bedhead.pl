#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: rename_methbed_head_hg18.pl
# Version: 0.1
# Last Updated: 5-6-2013
#
# This script will change the first line of a bed file to be something else
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input prefix
    2) NewID
" unless @ARGV == 2;

my $inprefixfile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $newid = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";




###################
# Another section #
###################

while(@Chr)
{
	
	my $firstline = <IN>;
    my $fline = chop($firstline);
    my @line = split ("PercMethylation", $fline);
    my @line2 = split ("PercentMethylation", $line[1]);
    print OUT $line[0] , $newid, $line2[0] , $newid , $line2[1], "\n";
	while (<IN>) { print OUT $_; }
	close IN;
	close OUT;
	

}

