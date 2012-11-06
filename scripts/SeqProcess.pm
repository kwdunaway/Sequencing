package SeqProcess;
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 10-30-2012
# Module Name: SeqProcess.pm
#
# This is a module with sequencing processing commands.
#
################################################################################################

###########################################################################
#                        Header and Path Modifier                         #
###########################################################################

sub add_path 
{
	my ($addtoPATH) = @_;

	`#!/bin/bash\n\n`;
	`PATH=\$PATH:$addtoPATH\n`;
	`export PATH\n\n`;
}

###########################################################################
#                 Combine and Filter Zipped Files                         #
#  Input: Raw file folder (only zipped files and the extension is .fq.gz) #
# Output: Filtered and Combined into one .fq file                         #
###########################################################################

sub filter_zip 
{
	my ($rawfqfolder) = @_;
	my $filtered_fastq = $rawfqfolder . "filtered.fq";

	`gunzip -c $rawfqfolder*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v \"^--\$\" >  $filtered_fastq\n\n`;

	return $filtered_fastq;
}
	
###########################################################################
#                               Run Bowtie                                #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) MM9 Path                                                     #
#         4) Filtered Fastq File                                          #
#                                                                         #
# Output: 1) Non-aligned Reads File                                       #
#         2) Aligned Preseparation File                                   #
###########################################################################

sub run_bowtie 
{
	my ($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq) = @_;
	
	`mkdir $ExperimentTopDir\n`;

	my $nonalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_NonAligned.fq";
	my $alignedpreseparationfile = $ExperimentTopDir . $BowtiePrefix . "_alignedpreseparation.txt";

	`bowtie -p 4 -M 1 -k 1 --chunkmbs 256 --strata --best --un $nonalignedreadsfile $mm9path $filtered_fastq $alignedpreseparationfile\n\n`;

	return ($nonalignedreadsfile, $alignedpreseparationfile);
}

###########################################################################
#                     Separate Repeats from Uniques                       #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) Aligned Preseparation File                                   #
#                                                                         #
# Output: 1) Unique Aligned Reads File                                    #
#         2) Repetitive Aligned Reads File                                #
###########################################################################

sub separate_repeats 
{
	my ($ExperimentTopDir, $BowtiePrefix, $alignedpreseparationfile) = @_;

	my $uniqalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Uniq.txt";
	my $repalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Repeat.txt";

	open(IN, "<$alignedpreseparationfile") or die "cannot open $alignedpreseparationfile infile";
	open(UNIQOUT, ">$uniqalignedreadsfile") or die "cannot open $uniqalignedreadsfile outfile";
	open(REPOUT, ">$repalignedreadsfile") or die "cannot open $repalignedreadsfile outfile";

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
			print UNIQOUT $_, "\n";
		}
	}


	close IN;
	close UNIQOUT;
	close REPOUT;

	return ($uniqalignedreadsfile, $repalignedreadsfile);
}

###########################################################################
#                     Eland Extended Format to BED                        #
#                                                                         #
#  This subroutine takes sequences in Eland Extended format (also known   #
#  as s_#_export.txt files) and produces multiple BED files which can     #
#  be used to analyze the data.                                           #
#                                                                         #
#  Input: 1) Input file name                                              #
#         2) Output file name                                             #
#         3) Read length                                                  #
#         4) Chromosome array number                                      #
#         5) Position array number                                        #
#         6) Strand array number                                          #
#         7) First Char                                                   #
#                                                                         #
# Output: Creates directory containing:                                   #
#           1) BED files for each chromosome                              #
#           2) Statistics on all reads                                    #
#           3) Gzipped non-unique and non-mappable reads files            #
#           4) Unknown reads file                                         #
###########################################################################

sub elandext_to_bed 
{
	# Input
	my ($infile, $outfile, $readlength, $chr, $pos, $strand, $firstchar) = @_;

	# Makes Output Directory
	if (! -d $outfile) 
	{ 
		`mkdir $outfile`; #creates dir if one doesn't exist
		if (! -d $outfile) { die "directory ($outfile) does not exist"} 
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
	my @ChrMapcount;
	for (my $n = 0; $n < 27; $n++)
	{
		$ChrMapcount[$n] = 0;
	}
	my $totalcount = 0; #Total Mapped Reads
	my $filename; # Outfile for each chromosome in bed
	my @chr_out; # Holds all the output files
	my $commandinput; # Holds Command Line Inputs

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
			if($array[$chr] =~ m/chr[0-9]/) # Chromosomes 1 to 23
			{
				for(my $n = 1; $n < 24; $n++)
				{	
					if($array[$chr] =~ m/(chr$n)$/)
					{ 	
						$chr_out[$n - 1]->print("chr$n \t" , $array[$pos] ,
								"\t" , $array[$pos]+$readlength, 									"\t", $outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
						$ChrMapcount[$n]++;
					}
				}
			}
			elsif($array[$chr] =~ m/chrX/) # Chromosome X
			{
				$chr_out[23]->print("chrX \t" , $array[$pos] , "\t" , 
							$array[$pos]+$readlength, "\t", 
							$outfile , "\t", "0", "\t" , 
							$array[$strand] , "\n");
				$ChrMapcount[24]++;
			}
			elsif($array[$chr] =~ m/chrY/) # Chromosome Y
			{
				$chr_out[24]->print("chrY \t" , $array[$pos] , "\t" , 
							$array[$pos]+$readlength, "\t", 
							$outfile , "\t", "0", "\t" , 
							$array[$strand] , "\n");
				$ChrMapcount[25]++;
			}
			elsif($array[$chr] =~ m/chrM/) # Chromosome of Mitochondria
			{
				$chr_out[25]->print("chrM \t" , $array[$pos] , "\t" , 
							$array[$pos]+$readlength, "\t", 
							$outfile , "\t", "0", "\t" , 
							$array[$strand] , "\n");
				$ChrMapcount[26]++;
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

}
# BED to WIG
# WIG to FPKMWIG
# Vis FPKMWIG



1;

