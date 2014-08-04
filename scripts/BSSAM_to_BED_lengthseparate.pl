#!/usr/bin/perl
use strict;
use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 7-28-2013
# Script Name: Remove_dup_readsv2.pl
#
# This script separates a SAM file into 2 bed files (based on given length of read).
#
# Arguments: See Below
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0.pl needs the following parameters:
    1) Input SAM file
    2) Read length separator (ex: 60)
    2) Above separate file bed file name
    3) Below and equal separate bed file name
" unless @ARGV == 4;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile IN infile";
my $sep_size = shift(@ARGV);
my $aboveoutfile = shift(@ARGV);
open(ABOVE, ">$aboveoutfile") or die "cannot open $aboveoutfile ABOVE outfile";
my $belowoutfile = shift(@ARGV);
open(BELOW, ">$belowoutfile") or die "cannot open $belowoutfile BELOW outfile";

while(<IN>){
	chomp;
	my @line = split("\t", $_);
	my $chr = $line[2];
	my $start = $line[3];
	my $end = $start + length($line[9]);
#	print $start , "\t" , $end , "\t", length($line[9]) ,"\n";
	if($end - $start > $sep_size){
		print ABOVE $chr , "\t" , $start , "\t", $end , "\n";
	}
	else{
		print BELOW $chr , "\t" , $start , "\t", $end , "\n";
	}
	
}
close IN;
close ABOVE;
close BELOW;

