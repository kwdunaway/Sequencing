#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 6-20-2012
# Script Name: splitfile.pl
#
# Splits a given file into subfiles with stated lines in each file.
#
# Arguments:
#    1) Input file
#    2) Output files prefix (everything before the number)
#    3) Output files suffix (everything after the number)
#    4) Number of lines per file
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "splitfile.pl needs the following parameters:
    1) Input file
    2) Output files prefix (everything before the number)
    3) Output files suffix (everything after the number)
    4) Number of lines per file
" unless @ARGV == 4;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outprefix = shift(@ARGV);
my $outsuffix = shift(@ARGV);
my $maxlines = shift(@ARGV);

my $filecounter = 1;
my $outfile = $outprefix . "00" . $filecounter . $outsuffix;
open(OUT, ">$outfile") or die "cannot open $outfile outfile";



#############
# Main Loop #
#############

my $linecounter = 0;
while(<IN>){
	if($linecounter % 4 ==0){
		my @check = split("", $_);
		#print $check[0] , "\n";
		if($check[0] eq "@"){}
		else {print "Skipping line: $linecounter in file $filecounter that looks like:\n$_\n"; next;}
		if($linecounter != 0) 
			{print "\n";}
		print OUT $_;
		my $line2 = <IN>;
		my $line3 = <IN>;
		my $line4 = <IN>;
		chop($line4);
		print $line2, $line3, $line4;
	}
	++$linecounter;
	if($linecounter > $maxlines){
		$linecounter = 1;
		close OUT;
		++$filecounter;
		if($filecounter < 10){
			$outfile = $outprefix . "00" . $filecounter . $outsuffix;
		}
		elsif($filecounter < 100){
			$outfile = $outprefix . "0" . $filecounter . $outsuffix;
		}
		else{
			$outfile = $outprefix . $filecounter . $outsuffix;
		}
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	}
}
		

