#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 6-10-2014
# Version: 1.1
#
# Takes a raw FASTQ file and filters out any reads found in the BS1 file by comparing
# The cell location.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following parameters:
    1) Input FASTQ file
    2) Input BS1 file
    3) Output FASTQ file
" unless @ARGV == 3;

my $FASTQinfile = shift(@ARGV);
open(FASTQ, "<$FASTQinfile") or die "cannot open $FASTQinfile infile";
my $BS1infile = shift(@ARGV);
open(BS1, "<$BS1infile") or die "cannot open $BS1infile infile";
my $outfile = shift(@ARGV);

my %FoundBS1;

print "Reading BS1 $BS1infile \n";
while(<BS1>){
	my @line = split("\t", $_);
	$FoundBS1{$line[0]} = 1;
	#HS1:273:C28V3ACXX:6:1107:7955:33170
	#HS1:273:C28V3ACXX:6:1101:1244:2243	
}
close BS1;

print "Filtering FASTQ $FASTQinfile \n";
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
while(<FASTQ>){
	my $line1 = $_;
	my $line2 = <FASTQ>;
	my $line3 = <FASTQ>;
	my $line4 = <FASTQ>;
	my @line = split(" ", $line1);
	my $ID =  substr $line[0], 1;
	if(defined $FoundBS1{$ID}) {}
	else {print OUT $line1, $line2, $line3, $line4;}	 
}
close FASTQ;
close OUT;
