#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Roy Chu
# Email: rgchu@ucdavis.edu
# Date: 1-8-2014
# Script Name: conversion_efficiency.pl
#
# This script calculates the conversion effiency of bisulfite sequencing. The input BED file 
# should be of mitochondrial DNA.
#
# Arguments:
#    1) Input BED file
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: conversion_efficiency.pl 
    1) Input BED file 
" unless @ARGV == 1;

my $inputname = shift(@ARGV);	# Input BED File
open(IN, "<$inputname") or die "conversion_efficiency.pl: Error: cannot open $inputname input BED file";
# Global Variables
my $totalreads = 0;
my $converted = 0;
my $unconverted = 0;
############################################
#           Reading Input BED File         #
############################################

#               Header of Input File    

my $header = <IN>;
if ($header =~ /track/)	# Header found
{	
	$header = <IN>; # Grab first line of data
}
else			# No header
{
}

# Process first line of data
chomp $header;
my @line = split("\t",$header);
my @percent = split("-", $line[3]);
$totalreads += $percent[1];			# Add to total reads
my $prerounded = $percent[0] * $percent[1]; 	# Multiply percentage unconverted by number of reads
my $rounded = int($prerounded + 0.5);		# Round to get unconverted reads
$unconverted += $rounded;			# Add to total unconverted reads
my $convrounded = $percent[1] - $rounded;	# Subtract to get converted reads
$converted += $convrounded;			# Add to total converted reads



############################################
#           Reading Input File             #
############################################

while(<IN>)
{ 
	chomp;
	my @line = split("\t",$_);
	my @percent = split("-", $line[3]);
	$totalreads += $percent[1];			# Add to total reads
	my $prerounded = $percent[0] * $percent[1]; 	# Multiply percentage unconverted by number of reads
	my $rounded = int($prerounded + 0.5);		# Round to get unconverted reads
	$unconverted += $rounded;			# Add to total unconverted reads
	my $convrounded = $percent[1] - $rounded;	# Subtract to get converted reads
	$converted += $convrounded;			# Add to total converted reads
}

my $efficiency = $converted/$totalreads;
print "Conversion Efficiency: ", $efficiency, "\n";
print "Total Reads: ", $totalreads, "\n";
print "Unconverted Reads: ", $unconverted, "\n";
print "Converted Reads: ", $converted, "\n";


close IN;


__END__

