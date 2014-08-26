#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway and Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 8-12-2014
#
# This script calculates average percent methylation of all CpG sites within each 
# line of a BED file. Multiple percent methylation folders with Percent Methylation
# bed files of each chromosome may be entered as inputs to be compared side by side.
#
# The user can set thresholds for each read. The minimum CpG site threshold 
# will place an "NA" for the read for that experiment if the specified amount
# of CpG sites found in that read is not met. The minimum read threshold will
# ignore CpG sites with reads pertaining to that site lower than the specified
# threshold. The minimum file threshold is useful when multiple folders are
# input and requires percent methylatiion data (not "NA") for a read from 
# at least the specified number of folders. If the file threshold is not met,
# the bed line is not printed to the output.
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: $0
    1) Output file
    2) Input BED or GTF File (needs to have a header line)
    3) Minimum CpG Site Threshold 
    4) Minimum Read Threshold
    5) Minimum File Threshold (Files without NA data)
    6+) Input Percent Methylation Folder Prefix (no chr)
" unless @ARGV > 5;

my $outputname = shift(@ARGV);	# Output with average percentage methylation per PMD
open(OUT, ">$outputname") or die "AvgMethv3.pl: Error: cannot open $outputname output file";
my $inputname = shift(@ARGV);	# Input BED or GTF File
open(IN, "<$inputname") or die "AvgMethv3.pl: Error: cannot open $inputname input BED file";
my $mincpg = shift(@ARGV);	# Avg % Meth = "NA" unless >= minimum CpG site threshold
# Special case (default): 0 as threshold will cause 0 CpG sites to give "NA"
my $minreads = shift(@ARGV);	# Threshold for reads found at a CpG site
my $minfiles = shift(@ARGV);	# Threshold for total data across all folders

# Global Variables
my %chrarray; 	   # Stores every field of lines in BED File and percentage methylation

############################################
#           Reading Input BED File         #
############################################
# Inputs all of the data in the BED file 

my $firstline = <IN>;		# Check for header
if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
	print "No header found, processing first line.\n";
	my @line = split("\t",$firstline);
	$chrarray{$line[0]}[0][0]++;	# Number of lines in this chromosome
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][0] = $line[0];	# Chromosome
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][1] = $line[1];	# Start
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][2] = $line[2];	# End
}
else { # If first line IS a header line
	print "Header found, skipping first line.\n";
}

# Process entire BED file
while(<IN>) { 
	chomp;
	my @line = split("\t",$_);
	if ($line[0] =~ /_/) {next;}
	$chrarray{$line[0]}[0][0]++;	# Number of lines in this chromosome
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][0] = $line[0];	# Chromosome
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][1] = $line[1];	# Start
	$chrarray{$line[0]}[$chrarray{$line[0]}[0][0]][2] = $line[2];	# End
}
close IN;
print "past loading CpG islands\n";

############################################
#           Reading Input Folders          #
############################################

# Run process and print for each chromosome
foreach my $key (sort keys %chrarray)
{
	my @outputtable;# Stores temporary output and prints at each iteration
#	my $count = 1;	# line count
	$outputtable[$count][0][2] = 0;		# Initialize number of CpG sites
	print "Loading $key...\n";
	for(my $i = 0; $i < $#ARGV+1; $i++)	# Run chromosome for each folder {
		# Open chromosome file
		my $filename = $ARGV[$i] . $key . ".bed";
		$count = 1;

		open (BED, "<$filename") or die "AvgMethv3.pl: Error: Couldn't open chromosome file $filename\n";
		print "Opening $filename \n";

		my $bedline = <BED>;  	# Check for header line
		if($bedline =~ /^chr/)	# Process line if data {
			my @bedlinearray = split("\t",$_);
			# If undefined, end reading of this file
			if(!defined $bedlinearray[0] || !defined $chrarray{$key}[$count][0]) {last;}
			# If diff chr, end reading of this file
			if($bedlinearray[0] ne $chrarray{$key}[$count][0])	{last;}
			# Go to next line if nothing found
			if($bedlinearray[2] > $chrarray{$key}[$count][2]) {
				if($outputtable[$count][$i][2] >= $mincpg && $outputtable[$count][$i][2] > 0){
					$outputtable[$count][0][3]++; # Number of files with info, used for minimum file threshold
				}
				$count++;
				$outputtable[$count][$i][2] = 0;
			} 
			else{
				# If within the range, then add to total meth
				if($bedlinearray[1] >= $chrarray{$key}[$count][1] && $bedlinearray[2] <= $chrarray{$key}[$count][2])
				{
					my @dash = split("-", $bedlinearray[3]);
					if ($dash[1] >= $minreads)
					{
						$outputtable[$count][$i][1] += $dash[0];	# Add percentage methylation
						$outputtable[$count][$i][2]++;	# Add to total sites
					}
				}
				# Otherwise, skip to starting location of next region to get in the range
			}			
		}


		while(<BED>)
		{
			my @bedlinearray = split("\t",$_);
			# If undefined, end reading of this file
			if(!defined $bedlinearray[0] || !defined $chrarray{$key}[$count][0])
			{
				last;
			}
			# If diff chr, end reading of this file
			if($bedlinearray[0] ne $chrarray{$key}[$count][0])	
			{
				last;
			}
			# Go to next line if nothing found
			if($bedlinearray[2] > $chrarray{$key}[$count][2]) 
			{
				# If read has data, add 1 to files with data
				if($outputtable[$count][$i][2] >= $mincpg && $outputtable[$count][$i][2] > 0)
				{
					$outputtable[$count][0][3]++; # Number of files with data for minimum file threshold
				}
				$count++; # next line
				$outputtable[$count][$i][2] = 0; #initialize
			} 
			else
			{
				# If within the range, then add to total meth
				if($bedlinearray[1] >= $chrarray{$key}[$count][1] && $bedlinearray[2] <= $chrarray{$key}[$count][2])
				{
					my @dash = split("-", $bedlinearray[3]);
					if ($dash[1] >= $minreads) # Read Threshold
					{
						$outputtable[$count][$i][1] += $dash[0];	# Add percentage methylation
						$outputtable[$count][$i][2]++;	# Add to total sites
					}
				}
				# Otherwise, skip to starting location of next region to get in the range
			}
		}
		close BED;

	}

	# Print data for this chromosome
	for(my $j = 1; $j <= $chrarray{$key}[0][0]; $j++)
	{
		# Check if files with data exist
		if(defined $outputtable[$j][0][3])
		{
			# Check minimum file threshold
			if($outputtable[$j][0][3] >= $minfiles)
			{
				# Print chromosome, start, end
				print OUT $chrarray{$key}[$j][0], "\t", $chrarray{$key}[$j][1], "\t", $chrarray{$key}[$j][2];

				# For each folder print data
				for(my $k = 0; $k < $#ARGV+1; $k++)
				{
					# Sites under CpG threshold -> "NA"
					if($outputtable[$j][$k][2] < $mincpg || $outputtable[$j][$k][2] == 0)	
					{
						print OUT "\tNA";
					}	
					# Sites Found -> Calculate average percent methylation	
					else				
					{
						$outputtable[$j][$k][0] = ($outputtable[$j][$k][1])/($outputtable[$j][$k][2]); # Calculate average
						my $rounded = sprintf("%.5f", $outputtable[$j][$k][0]);	
						print OUT "\t", $rounded;
					}
				}
				print OUT "\n";
			}
		}
	}
}

close OUT;
