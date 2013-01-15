#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: bed_to_userdefGTF.pl
# Version: 1.0
# Last Updated: 12-11-2012
#
# This script convert a bed file to a userdef GTF file that is able to be used 
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Bed infile name
    2) GTF outfile name
    3) Track description
" unless @ARGV == 3;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
my $filedesc = shift(@ARGV);



###################
# Another section #
###################

<IN>;  # take of header line
my $inlinenum = 1;
while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    if($inlinenum > 1) 
    {print OUT "\n";}
    my $blockname = "Chunk_";
    if($inlinenum < 10){$blockname = $blockname . "000" . $inlinenum;}
    elsif($inlinenum < 100){$blockname = $blockname . "00" . $inlinenum;}
    elsif($inlinenum < 1000){$blockname = $blockname . "0" . $inlinenum;}
	else{$blockname = $blockname . $inlinenum;}
    print OUT $line[0], "\t" , $filedesc , "\t" , $blockname , "\t", $line[1], "\t", $line[2], "\t.\t+";
    $inlinenum++;
}

