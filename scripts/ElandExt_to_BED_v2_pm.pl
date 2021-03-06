#!/usr/bin/perl -w
use strict; use warnings;

# ElandExt_to_BEDv2.pl          6/6/12 version
#
# This program takes sequences in Eland Extended format (also known as s_#_export.txt files) and 
# produces multiple files BED which can be used to analyze the data.  The bed files go out to the 
# 6th column (strand notation) but have null values for columns 4 (name) and 5 (score).  This
# program will also give statistics and separate out nonuniquely mapped reads.
#
# This is a little more flexible than other ElandExt_to_BED versions because all you need for 
# chromosomal notation is chr##  (ex: chr1 or chr13) somewhere in the chrom notation.  Also,
# you can input read sequence length (which is better than the program figuring it out for speed
# purposes.  Finally, you must input which array field the chromosome, position, and strand are
# in just in case they are different from export file to export file (10, 11, 12 for Berkeley files).
#
# Arguments:
#   1) Input file name
#   2) Output file name
#   3) Read length
#   4) Chromosome array number
#   5) Position array number
#   6) Strand array number
#   7) First Char



##################################################
# Command Line Error Checking and I/O Initiation #
##################################################

die "useage: ElandExt_to_BEDv2.pl 
<input file name> 
<output filename>
<read length>
<Chromosome array number>
<Position array number>
<Strand array number>
<First Char>" unless @ARGV == 9;
my $infile = shift(@ARGV);
my $ExperimentTopDir = shift(@ARGV); 
my $FilePrefix = shift(@ARGV); 
my $basereadlength = shift(@ARGV); 
my $finalreadlength = shift(@ARGV); 
my $chr = shift(@ARGV); 
my $pos = shift(@ARGV); 
my $strand = shift(@ARGV); 
my $MaxDupReads = shift(@ARGV); 


 	print "\nBeginning conversion of Eland Extended format to BED format\n";

	my $outdir = $ExperimentTopDir . "/" . $FilePrefix . "_bed";
	my $minusstrandlength = $finalreadlength - $basereadlength;
	my %Count;
	my %Files;
	my $totalcount = 0;

	# Makes Output Directory
	print "Making $outdir directory\n";
	if (! -d $outdir) 
	{ 
		`mkdir $outdir`; #creates dir if one doesn't exist
		if (! -d $outdir) { die "directory ($outdir) does not exist"} 
	}

	open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read
	print "Processing data into each chromosome output file in BED format\n";
	while (<IN>)
	{
		chomp;
		my @array = split("\t", $_); # Splitting data into array
		my $chrom = $array[$chr];
		my $readstrand = $array[$strand];
		$totalcount++; # Add to total mapped reads

		#If outputfile has not been for the chromosome, make it
		if(!$Count{$totalcount}){
			$Count{$chrom} = 0;
			my $filename = $outdir . "/" . $FilePrefix . "_" . $chrom . ".bed";
			open($Files{$chrom}, ">$filename") or die "cannot open $filename outfile";	
		}
		my $printfile = $Files{$chrom};
		$Count{$chrom}++;

		if($readstrand == "+"){
			print {$Files{$chrom}} $chrom , "\t" , $array[$pos] , "\t" , $array[$pos] + $finalreadlength, "\t", $FilePrefix , "\t", "0", "\t" , $readstrand , "\n";
		}
		elsif($readstrand == "-"){
			print {$Files{$chrom}} $chrom , "\t" , $array[$pos] - $minusstrandlength, "\t" , $array[$pos] + $basereadlength, "\t", $FilePrefix , "\t", "0", "\t" , $readstrand , "\n";
		}
		else {die "Strand is not + nor -, it is $readstrand";}
	}
	close IN;
	print "Finished outputting data to each BED file\n";


	############################################################################
	#                  Printing statistics to Stats Outfile                    #
	############################################################################

	my $filename = $ExperimentTopDir . "/" . "Stats_" . $FilePrefix . ".txt";
	print "Printing statistics to ", $filename, "\n";
	open(STATS, ">$filename") or die "cannot open $filename outfile";	
	foreach my $key (sort keys %Files) {
		print STATS "Number of reads mapped " , $key , " is:\t" , $Count{$key} , "\n";
		close $Files{$key};
		my $bedfile = $outdir . "/" . $FilePrefix . "_" . $key . ".bed";
		sort_bed($bedfile);
		eliminate_bed_dups($bedfile, $MaxDupReads);
	}
	print STATS "\nTotal Number of Total mapped reads is:\t", $totalcount, "\n";
	close STATS;
	print "Finished conversion of Eland Extended to BED format\n";

sub sort_bed
{
	my ($bedfile) = @_;
	my $temp = $bedfile . "_sorted.bed";
	`sort -n +1 -2 $bedfile > $temp`;
	`rm $bedfile`;
	`mv $temp $bedfile`;
}

