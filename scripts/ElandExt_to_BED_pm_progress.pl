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
<First Char>" unless @ARGV == 7;
my $infile = shift(@ARGV);
my $outfile = shift(@ARGV); 
my $readlength = shift(@ARGV); 
my $chr = shift(@ARGV); 
my $pos = shift(@ARGV); 
my $strand = shift(@ARGV); 
my $firstchar = shift(@ARGV); 
if (! -d $outfile)
{ 
     `mkdir $outfile`; #creates dir if one doesn't exist
     if (! -d $outfile) { die "directory ($outfile) does not exist"} #dir still does not exist
}



	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read

	my @array; # Contains the data to be parsed
	my $QCcount = 0; # Low Quality Reads (placed in non-mappable reads)
	my $NMcount = 0; # Non-mappable Reads
	my $NonUniqcount = 0; # Non-unique Reads
	my $Unknowncount = 0; # Unknown Reads
	my @ChrMapcount; # Unique Mapped Reads by Chromosome
	for (my $n = 0; $n < 27; $n++)
	{
		$ChrMapcount[$n] = 0;
	}
	my $totalcount = 0; #Total Mapped Reads
	my $filename; # Outfile for each chromosome in bed
	my @chr_out; # Holds all the output files
	my $commandinput; # Holds Command Line Inputs

	# Hash made to include X, Y, and M
	my %chromosomes = (); # Contains chromosomes(keys = call number, values = chromosome number)

	for(my $n = 1; $n < 24; $n++)
	{
		$chromosomes{$n} = $n;
	}
	$chromosomes{24} = "X";
	$chromosomes{25} = "Y";
	$chromosomes{26} = "M";

	############################################################################
	#       Creates BED file for each Chromosome, Puts Output files in array   #
	#                                                                          #
	#   Chromosome 1 corresponds to $chr_out[0], Chr2 -> $chr_out[1], etc.     #
	#                    ChrX = [23], ChrY = [24], ChrMitochondria = [25]      #
	#   *Fastq Files* -> Non-mappable = [26], Non-unique = [27]                #
	#   *Text Files*  -> Stats = [28], Unknown = [29]                          #
	############################################################################

	for(my $n = 0; $n < 30; $n++)
	{
		if($n > 22)
		{
			if($n == 23)
			{
				$filename = $outfile . "/" . $outfile . "_chrX.bed";
			}
			if($n == 24)
			{
				$filename = $outfile . "/" . $outfile . "_chrY.bed";
			}
			if($n == 25)
			{
				$filename = $outfile . "/" . $outfile . "_chrM.bed";
			}
			if($n == 26)
			{
				$filename = $outfile . "/" . $outfile . "_NM.fq"; 
			}
			if($n == 27)
			{
				$filename = $outfile . "/" . $outfile . "_NonUnique.fq";
			}
			if($n == 28)
			{
				$filename = $outfile . "/" . "Stats_" . $outfile . ".txt";
			}
			if($n == 29)
			{
				$filename = $outfile . "/" . "Unknown_" . $outfile . ".txt";
			}
		}
		else
		{
			$filename = $outfile . "/" . $outfile . "_chr" . ($n + 1) . ".bed";
		}
		open($chr_out[$n], ">$filename") or die "cannot open $filename outfile";
	}


	############################################################################
	#           Processing data into each chromosome output file               #
	############################################################################

	while (<IN>)
	{
		chomp;
		@array = split("\t", $_); # Splitting data into array
		$totalcount++; # Total mapped reads
		# Non-mappable reads
		if ($array[$chr] eq "QC") # Low quality reads (placed in non-mappable reads)
		{ 
			$chr_out[26]->print("$_" , "\n"); 
			$QCcount++;
		}
		elsif ($array[$chr] eq "NM") # Non-mappable reads
		{ 	
			$chr_out[26]->print("$_" , "\n");  
			$NMcount++;
		}
		# Non-unique reads
		elsif ($array[$chr] =~ m/:/)
		{ 	
			$chr_out[27]->print("$_" , "\n"); 
			$NonUniqcount++;
		}
		# Processing unique reads to chromosomes
		elsif($array[$chr] =~ m/chr/)
		{
			foreach my $chromosome(keys %chromosomes) # Check all chromosomes
			{	
				my $chr_value = $chromosomes{$chromosome};
				if($array[$chr] =~ m/(chr$chr_value)$/)
				{ 	
					$chr_out[$chromosome-1]->print("chr", $chr_value, "\t" , 
								$array[$pos] , "\t" , 
								$array[$pos]+$readlength, 									"\t", $outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[$chromosome]++;
				}
			}
			$ChrMapcount[0]++; # Add to total number of uniquely mapped reads
		}
		else # Unknown data
		{
			$chr_out[29]->print("$_\n");
			$Unknowncount = $Unknowncount + 1;
  		}
	}

	############################################################################
	#                  Printing statistics to Stats Outfile                    #
	############################################################################

	my $weirdcount = $totalcount - $QCcount - $NMcount - $NonUniqcount - $ChrMapcount[0];

	for (my $n = 1; $n < 24; $n++)
	{
		if ($ChrMapcount[$n] > 0)
		{
			$chr_out[28]->print("Number of reads mapped to Chromosome $n is:\t" , 
								$ChrMapcount[$n], "\n");
		}
	}
	$chr_out[28]->print("Number of reads mapped to Chromosome X is:\t", $ChrMapcount[24], "\n");
	$chr_out[28]->print("Number of reads mapped to Chromosome Y is:\t", $ChrMapcount[25], "\n");
	$chr_out[28]->print("Number of reads mapped to Mitochondria is:\t", $ChrMapcount[26], "\n\n");
	$chr_out[28]->print("Number of Uniquely mapped reads is:\t", $ChrMapcount[0], "\n");
	$chr_out[28]->print("Number of NonUniquely mapped reads is:\t", $NonUniqcount, "\n");
	$chr_out[28]->print("Number of QC (Low Quality reads) is:\t", $QCcount, "\n");
	$chr_out[28]->print("Number of NM (nonmappable reads) is:\t", $NMcount, "\n");
	$chr_out[28]->print("Number of unknown results (printed in separate file) is:\t", 
								$Unknowncount, "\n");
	$chr_out[28]->print("Number of weird reads (should be same as previous line) is:\t", 
								$weirdcount, "\n");
	$chr_out[28]->print("Number of Total mapped reads is:\t", $totalcount, "\n");

	############################################################################
	#      Closing Files and Gzip Non-mappable and Non-unique Reads Files      #
	############################################################################

	# Closing files
	for(my $n = 0; $n < 30; $n++)
	{
		close($chr_out[$n]);
	}
	close(IN);

	# Gzip Non-mappable and non-unique reads outfiles
	$commandinput = "gzip " . $outfile . "/" . $outfile . "_NM.fq";
	`$commandinput`;
	$commandinput = "gzip " . $outfile . "/" . $outfile . "_NonUnique.fq";
	`$commandinput`;

	############################################################################
	#                            Sort All Bed Files                            #
	############################################################################

	my $bedfile;

	foreach my $chromosome(keys %chromosomes) # Check all chromosomes
	{
		my $chr_value = $chromosomes{$chromosome};
		$bedfile = $outfile . "/" . $outfile . "_chr" . $chr_value . ".bed";
		if ($ChrMapcount[$chromosome] == 0) # Check if empty
		{
			`rm $bedfile`; # Delete if empty
		}
		else
		{
			sort_bed($bedfile);
		}
	}

sub sort_bed
{
	my ($bedfile) = @_;
	my $temp = $bedfile . "_sorted.bed";
	`sort -n +1 -2 $bedfile > $temp`;
	`rm $bedfile`;
	`mv $temp $bedfile`;
}

