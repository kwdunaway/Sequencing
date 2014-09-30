#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 9-15-2014
#
# This script will take the input prefix of a percent methylation file and output
# a "typical" Methylation score file as described by methylSig. The file has the
# following format:
#
#  	chrBase   chr     base strand coverage  freqC freqT
#  	chr21.43008527 chr21 43008527      F	32 100.00  0.00
#
# Since the information for forward and reverse has been combined for CpG islands,
# the file will only have F strand.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input prefix (leave off ##.bed)
    2) Output prefix
" unless @ARGV == 2;

my $inprefix = shift(@ARGV);
my $outprefix = shift(@ARGV);
#open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";

my @chrom;
for (my $n = 1; $n < 23; $n++){
	push(@chrom, $n);
}
push(@chrom, "X");
push(@chrom, "Y");

###################
# Another section #
###################

while(@chrom)
{
	my $chr = shift(@chrom);
	my $infile = $inprefix . $chr . ".bed";
	open(IN, "<$infile") or die "cannot open $infile infile";
	print "Starting $infile\n";
	my $outfile = $outprefix . "chr" . $chr . ".mscore";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";
	
	#remove header
	<IN>;
	#print new header
	print OUT "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";

	while (<IN>)
	{
	    chomp;
	    my @line = split ("\t", $_);
	    my $name = $line[0] . "\." . $line[1];
	    my @methinfo = split ("-", $line[3]);
	    my $coverage = $methinfo[1];
	    my $freqC = sprintf("%.2f", 100*$methinfo[0]);
	    my $freqT = sprintf("%.2f", 100*(1-$methinfo[0]));
	    print OUT $name,"\t",$line[0],"\t",$line[1],"\tF\t",$coverage,"\t",$freqC,"\t",$freqT,"\n";
	}
	close IN;
	close OUT;
}
