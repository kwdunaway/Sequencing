#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 9-29-2014
#
# Takes SAM file(s) as input, removes duplicate reads (keeps the longest) and combines
# them all into a single file. This file must then be sorted using samtools sort.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Outfile
    2+) Input Sorted SAM file
" unless @ARGV > 1;

my $outfile = shift(@ARGV);
#open(OUT, ">$outfile") or die "cannot open $outfile outfile";
open(OUT, ">tmp.sam") or die "cannot open tmp.sam outfile";
my @Infiles = shift(@ARGV);

#Columns for formatting SAM files
my $chrc = 2;
my $startc = 3;
my $methc = 14;
my $strandc = 11;



#############
# Main Loop #
#############

while (@Infiles){
	my $infile = shift(@Infiles);
	open(IN, "<$infile") or die "cannot open $infile infile";

	# Previous read variables
	my $prevline = <IN>;
	chop($prevline);
	my @line = split ("\t", $prevline);
	my $chrom = $line[$chrc];
	my $start = $line[$startc];
	my $methstring = substr $line[$methc], 5;
	my $strand = substr $line[$strandc], 5,1;
	my $prevmethstring = $methstring;
	my $prevstart = $start;
	my $prevstrand = $strand;
	
#	die "$line[$methc] \n ";

	while (<IN>)
	{
	    chomp;
	    @line = split ("\t", $_);
		$chrom = $line[$chrc];
		$start = $line[$startc];
		$methstring = substr $line[$methc], 5;
		$strand = substr $line[$strandc], 5,1;

		# If duplicate line, take longest read and then skip
		if($prevstart == $start && $prevstrand eq $strand) {
			if(length($methstring) > length($prevmethstring)){
				$prevline = $_;
				$prevmethstring = $methstring;
			}
		}
		else{
			print OUT $prevline , "\n";
			$prevline = $_;
			$prevmethstring = $methstring;
			$prevstart = $start;
			$prevstrand = $strand;
		}
	}
	close IN;
}
close OUT;
my $commandline = "sort --key=3,3 --key=4,4n tmp.sam > " . $outfile;
`$commandline`;
`rm tmp.sam`;
