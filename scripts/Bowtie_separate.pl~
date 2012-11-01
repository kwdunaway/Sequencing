#!/usr/bin/perl -w
use strict; use warnings;

# Filter.pl 4/15/11 version
#
# This is a simple filter program that will take a bowtie outupt and
# separate it into two files based on the 7th column ($line[6]).  It
# will be used to separate repetative aligned reads from non repetative
# ones.
#
# Arguments:
#   1) Infile name
#   2) Non Rep Outfile name
#   3) Repetative Outfile name



##################################################
# Command Line Error Checking and I/O Initiation #
##################################################

die "useage: Bowtie_separate.pl 
   1) Infile name
   2) Non Rep Outfile name
   3) Repetative Outfile name
" unless @ARGV == 3;
my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $nonoutfile = shift(@ARGV);
open(NONOUT, ">$nonoutfile") or die "cannot open $nonoutfile  outfile";
my $repoutfile = shift(@ARGV);
open(REPOUT, ">$repoutfile") or die "cannot open $repoutfile  outfile";

print "\n Filtering $infile to: \n Nonrep Outfile named $nonoutfile \n and Rep Outfile named $repoutfile \n";

 
while (<IN>)
{
    chomp;
    my @line = split ("\t", $_);
    if($line[6] > 0)
    {
	print REPOUT $_, "\n";
    }
    else
    {
	print NONOUT $_, "\n";
    }
}


close IN;
close NONOUT;
close REPOUT;
