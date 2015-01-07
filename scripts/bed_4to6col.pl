#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 11-25-2014
#
# This script takes in a 4 column bed file and adds 2 more columns. The added columns 
# will be for name and strand. The 4th column in the input file will be moved
# to the score column (because we are assuming they are paired read lengths which
# are between 0 and 1000.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Input bed file (4 columns) (can read zipped files)
    2) Output bed file (6 columns) (can print zipped or unzipped files)
" unless @ARGV == 2;

my $infile = shift(@ARGV);
if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
else{open(IN, "<$infile") or die "cannot open $infile infile";}
my $outfile = shift(@ARGV);
my $outfile_no_name = $outfile;
if ($outfile_no_name =~ /\.gz$/) {
	$outfile_no_name = substr($outfile_no_name, 0 , -3);
}
open(OUT, ">$outfile_no_name") or die "cannot open $outfile_no_name outfile";


###################
# Another section #
###################

my $read_number = 0;

while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t" , "READ_" , $read_number, "\t", $line[3], "\t+\n";
    $read_number++;
}
close IN;
close OUT;

if ($outfile =~ /\.gz$/) {
	my $commandline = "gzip " . $outfile_no_name;
	print "Zipping $outfile_no_name\n";
	`$commandline`;
}
