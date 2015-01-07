#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 12-8-2014
#
# This script does multiple processes to create two fastq files from an already quality
# filtered and combined fastq file.
#
# 1) It separates fastq files for those with and without adapter sequence. If it doesn't
# contain adapter contamination, those will be put into a file unprocessed.
#
# 2.1) If there is adapter contamination, it trims adapter sequence from a fastq file.
# Currently, it only takes the first 10 bases of the adapter sequence, searches for it, 
# then trims any read that has a full match for the adapter.
#
# 2.2) This script also "chews" back (removes) the last X (user defined) bases at the 3' 
# end. This was included because we found many bases that close to the end to have
# a unmethylated skew. 
#
# 2.3) The minimum read length checks after trimming and chewing occurs. If the read does
# not meet this length, it gets removed from the dataset.
#
# 3) There is also a filtering of reads with a quality score containing "####". This
# was done to quickly remove low quality reads that usually are a result of 
# adapter dimers.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "Usage: $0 needs the following parameters:
    1) Input FASTQ file to be split (can read uncompressed or gzipped files)
    2) Output FASTQ file name (no orig adapter contam) (Uncompressed unless .gz is at the end of the file name)
    3) Output FASTQ file name (trimmed adapter) (Uncompressed unless .gz is at the end of the file name)
    4) Minimum read length after adapter trimming (ex: 35)
    5) Chew back length (ex: 10)
" unless @ARGV == 5;

my $infile = shift(@ARGV);
if ($infile =~ /\.gz$/) {open(IN, "gunzip -c $infile |") or die "can't open pipe to $infile";}
else{open(IN, "<$infile") or die "cannot open $infile infile";}
#open(IN, "<$infile") or die "cannot open $infile infile";
my $outfile_no = shift(@ARGV);
my $outfile_no_name = $outfile_no;
if ($outfile_no_name =~ /\.gz$/) {
	$outfile_no_name = substr($outfile_no_name, 0 , -3);
}
open(OUTN, ">$outfile_no_name") or die "cannot open $outfile_no_name outfile";
my $outfile_with = shift(@ARGV);
my $outfile_with_name = $outfile_with;
if ($outfile_with_name =~ /\.gz$/) {
	$outfile_with_name = substr($outfile_with_name, 0 , -3);
}
open(OUTA, ">$outfile_with_name") or die "cannot open $outfile_with_name outfile";
my $minreadlength = shift(@ARGV);
my $trimlength = shift(@ARGV);

#Solexa_reverse_contam
#AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAA
#my $adapter = "AGATCGGAAG";
my $fulladapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAA";

#my %Sequences_After_Trimming = Print_Trimmed_FASTQ_file($infile, $outfile, $fulladapter, $minreadlength, $trimlength);
my $adapter_seq = substr($fulladapter, 0 , 10);
my $counter = 0;



####################
# Main Infile Loop #
####################

while (<IN>) {
	# Take in 4 lines, aka a single read
	chomp;
    my $ID = $_;
  	my $seq = <IN>;
 	chop($seq);
	my $third = <IN>;
	$third = "+";
	my $quality = <IN>;
    chop($quality);
    
	# Processing message
	$counter++;
	if($counter % 1000000 == 0) {print "Finished procesing read:\t" , $counter , "\n";}

	# Quality control
	if($quality =~ m/####/) {next;}
    if($seq =~ m/$adapter_seq/) {
		# Trim sequence (because of adapter contamination)
    	my @trimmedseq = split($adapter_seq, $seq);
		$seq = $trimmedseq[0];
		my $seqlength = length($seq) - $trimlength;
		$seq = substr($seq, 0, $seqlength);
	    $quality = substr($quality, 0, $seqlength);
		if($seqlength >= $minreadlength) {
	    	print OUTA $ID , "\n" , $seq , "\n" , $third , "\n" , $quality , "\n";
	    }
	}
	else{
		# No adapter contaminated sequences
	    print OUTN $ID , "\n" , $seq , "\n" , $third , "\n" , $quality , "\n";
	}
}
close IN;
close OUTA;
close OUTN;


################################
# Outfile zipping if necessary #
################################

if ($outfile_no =~ /\.gz$/) {
	my $commandline = "gzip " . $outfile_no_name;
	print "Zipping $outfile_no_name\n";
	`$commandline`;
}
if ($outfile_with =~ /\.gz$/) {
	my $commandline = "gzip " . $outfile_with_name;
	print "Zipping $outfile_with_name\n";
	`$commandline`;
}
