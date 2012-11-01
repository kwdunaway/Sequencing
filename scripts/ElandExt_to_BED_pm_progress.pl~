#!/usr/bin/perl -w
use strict; use warnings;

# ElandExt_to_BEDv2.pl          6/6/12 version
#
# This program takes sequences in Eland Extended format (also known as s_#_export.txt files) and 
# produces multiple files BED which can be used analyze the data.  The bed files go out to the 
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

	my @array;
	my $QCcount = 0;
	my $NMcount = 0;
	my $NonUniqcount = 0;
	my $Unknowncount = 0;
	my @ChrMapcount;
	for (my $n = 0; $n < 27; $n++)
	{
		$ChrMapcount[$n] = 0;
	}
	my $totalcount = 0;
	my $filename; # Outfile for each chromosome in bed
	my @chr_out; # Holds all the output files

	############################################################################
	#       Creates BED file for each Chromosome, Puts Output files in array   #
	#                                                                          #
	#   Chromosome 1 corresponds to $chr_out[0], Chr2 -> $chr_out[1], etc.     #
	#                    ChrX = [23], ChrY = [24], ChrMitochondria = [25]      #
	#   *Fastq Files* -> Non-mappable = [26], Non-unique = [27]                  #
	#   *Text Files*  -> Stats = [28], Unknown = [29]                          #
	############################################################################

	for(my $n = 0; $n < 30; $n++)
	{
		if($n > 22)
		{
			if($n == 23)
			{
				$filename = $outfile . "/" . $outfile . "_chrX" . ".bed";
			}
			if($n == 24)
			{
				$filename = $outfile . "/" . $outfile . "_chrY" . ".bed";
			}
			if($n == 25)
			{
				$filename = $outfile . "/" . $outfile . "_chrM" . ".bed";
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
		@array = split("\t", $_);
		$totalcount++;
		# Non-mappable Reads
		if ($array[$chr] eq "QC")
		{ 
			$chr_out[26]->print("$_" , "\n"); 
			$QCcount++;
		}
		elsif ($array[$chr] eq "NM")
		{ 	
			$chr_out[26]->print("$_" , "\n");  
			$NMcount++;
		}
		# Non-unique Reads
		elsif ($array[$chr] =~ m/:/)
		{ 	
			$chr_out[27]->print("$_" , "\n"); 
			$NonUniqcount++;
		}
		else 
		{
			for(my $n = 1; $n < 24; $n++)
			{
				if($array[$chr] =~ m/chr$n/)
				{ 	
					$chr_out[$n - 1]->print("chr$n \t" , $array[$pos] , "\t" , 
								$array[$pos]+$readlength, "\t", 
								$outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[$n]++;
				}
				elsif($array[$chr] =~ m/chrX/)
				{
					$chr_out[23]->print("chrX \t" , $array[$pos] , "\t" , 
								$array[$pos]+$readlength, "\t", 
								$outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[24]++;
				}
				elsif($array[$chr] =~ m/chrY/)
				{
					$chr_out[24]->print("chrY \t" , $array[$pos] , "\t" , 
								$array[$pos]+$readlength, "\t", 
								$outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[25]++;
				}
				elsif($array[$chr] =~ m/chrM/)
				{
					$chr_out[25]->print("chrM \t" , $array[$pos] , "\t" , 
								$array[$pos]+$readlength, "\t", 
								$outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[26]++;
				}
				$ChrMapcount[0]++;
			}
		}
	}

	for(my $n = 0; $n < 30; $n++)
	{
		close($chr_out[$n]);
	}
	close(IN);


