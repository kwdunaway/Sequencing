#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: Line1_analysis.pl
# Version: 0.1
# Last Updated: April 222, 2013
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the Line1 pattern.  Then, it puts those reads in a separate file.  Also, the
# script quantifies methylation of these sequences across three potential 
# methylation sites.
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input fastq file
    2) Output fastq file with Line1 reads
    3) Stats file
" unless @ARGV == 3;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";
my $stats_filename = shift(@ARGV);
open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";


###################
# Another section #
###################

my $site1 = 0;
my $site2 = 0;
my $site3 = 0;
my $LINE1 = 0;

while (<IN>)
{
    chomp;
    if($_ =~ /TTYGTGGTG[CT]GT[CT]GTTTTTT/){
    	$LINE1++;
    	print OUT $_ , "\n";
    }
}


__END__
TTYGT
GGTGY
GTYGT
TTTTT
A(A(t?)G)TY
GGTTT

