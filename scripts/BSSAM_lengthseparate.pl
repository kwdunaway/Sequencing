#!/usr/bin/perl
use strict;
use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 9-9-2014
#
# This script separates a SAM file into 2 sam files (based on given length of read).
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
    2) Above separate file sam file name
    3) Below and equal separate sam file name
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
	if(length($line[9]) > $sep_size){
		print ABOVE $_, "\n";
	}
	else{
		print BELOW $_, "\n";
	}
	
}
close IN;
close ABOVE;
close BELOW;

