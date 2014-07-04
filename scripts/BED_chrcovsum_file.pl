#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 6-20-2014
# Version: 1.0
#
# Takes a bed file and sums the coverage (in bases) of each chromosome. Prints to screen.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "BED_chrcovsum needs the following parameters:
    1) Input BED file to be summed
" unless @ARGV == 1;

my $infile = shift(@ARGV);

my %Results = BED_chrcovsum_file($infile);
foreach my $chr (sort(keys %Results)){
	print $chr , "\t", $Results{$chr} , "\n";
}



###############
# Subroutines #
###############

sub BED_chrcovsum_file{
	my ($BEDinfile) = @_;
	my %BEDhash;
	my %Resultshash;
	open(BEDIN, "<$BEDinfile") or die "cannot open $BEDinfile BEDinfile";
	while(<BEDIN>){
		my @line = split("\t", $_);
		my $chr = $line[0];
		if (substr($chr,0,3) ne "chr") {next;}
		my $start = $line[1];
		my $stop = $line[2];
		$BEDhash{$chr}{$start} = $stop;		
	}
	foreach my $chr (sort(keys %BEDhash)){
		my $runningtotal = 0;
		foreach my $start (keys %{$BEDhash{$chr}}) {
			$runningtotal = $runningtotal + $BEDhash{$chr}{$start} - $start;
		}
		$Resultshash{$chr} = $runningtotal;
	}
	return(%Resultshash)
}
