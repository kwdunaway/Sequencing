#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 11/3/2014
#
# This script filters 3 of same character kmers out
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Input file name
    2) output file name
" unless @ARGV == 2;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";


###################
# Another section #
###################

while (<IN>)
{
    chomp;
    if($_ =~ /TTT/){next;}
    elsif($_ =~ /AAA/){next;}
    elsif($_ =~ /GGG/){next;}
    elsif($_ =~ /CCC/){next;}
    print OUT $_ , "\n";
}
close IN;
close OUT;
