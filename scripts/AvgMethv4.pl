#!/usr/bin/perl 
use strict; use warnings;

################################################################################
# Author: Keith Dunaway and Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 8-14-2014
#
# This script calculates average percent methylation of all CpG sites in each 
# line of a BED file. Multiple percent methylation folders with sorted
# percent methylation bed files of each chromosome may be entered as inputs to
# be compared side by side.
#
# The user can set thresholds for each read. The minimum CpG site threshold 
# will place an "NA" for the read for that experiment if the specified amount
# of CpG sites found in that read is not met. The minimum read threshold will
# ignore CpG sites with reads pertaining to that site lower than the specified
# threshold. The minimum file threshold is useful when multiple folders are
# input and requires percent methylation data (not "NA") for a read from 
# at least the specified number of folders. If the file threshold is not met,
# the bed line is not printed to the output.
#
# Arguments:
#    <see below>
#
################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: $0
    1) Output file
    2) Input BED or GTF File (needs to have a header line)
    3) Minimum CpG Site Threshold 
    4) Minimum Read Threshold
    5) Minimum File Threshold (Files without NA data)
    6,8+) Input Percent Methylation Folder Prefix (no chr)
    7,9+) Input Sample Name (for header of output file)
" unless @ARGV > 5;

my $outputname = shift(@ARGV);	# Output with average percentage methylation per PMD
open(OUT, ">$outputname") or die "AvgMethv3.pl: Error: cannot open $outputname output file";
my $inputname = shift(@ARGV);	# Input BED or GTF File
open(BED, "<$inputname") or die "AvgMethv3.pl: Error: cannot open $inputname input BED file";
my $mincpg = shift(@ARGV);	# Avg % Meth = "NA" unless >= minimum CpG site threshold
# Special case (default): 0 as threshold will cause 0 CpG sites to give "NA"
my $minreads = shift(@ARGV);	# Threshold for reads found at a CpG site
my $minfiles = shift(@ARGV);	# Threshold for total data across all folders
my @filenames;			# Input Percent Methylation Folder Prefix (Ex: /home/user/Permeth/Permeth_)
my @headernames;		# Input Percent Methylation Header for Column (Ex: Sample1)
while(@ARGV){
	push(@filenames,shift(@ARGV));
	push(@headernames,shift(@ARGV));
}

# Global Variables
my %bed_hash;	# Stores every field of lines in BED File and percentage methylation
my %bed_array;	# Array version of the hash

##########################################################################
#                       Reading Input BED File                           #
# First, the BED file is loaded onto a hash by first checking for        #
# a header and then splitting every line by chromosome, start, and stop. #
# The hash is then sorted into an array to be used in later steps for    #
# faster data comparison.                                                #
##########################################################################

my $firstline = <BED>;		# Check for header
if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
	print "No header found, processing first line.\n";
	chomp($firstline);
	my @line = split("\t",$firstline);
	$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Push line to hash
}
else { # If first line IS a header line
	print "Header found, skipping first line!\n";
}

# Process entire BED file
while(<BED>)
{ 
	chomp;
	my @line = split("\t",$_);
	if ($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
	# Check if duplicate
	if(!defined $bed_hash{$line[0]}{$line[1]}{$line[2]}) {
		$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Push line to hash
	}
	# Otherwise duplicate
	else {
		print "Warning: duplicate found, skipping: $_\n";
	}
}

# Fill in array with hash information, this is done to sort the information
foreach my $chr (keys %bed_hash){
	# Initialize line count
	$bed_array{$chr}[0][2]=0;
	# For each position, sort numerically
	foreach my $start (sort {$a<=>$b} keys %{$bed_hash{$chr}}){
		foreach my $end (sort {$a<=>$b} keys %{$bed_hash{$chr}{$start}}){
			$bed_array{$chr}[$bed_array{$chr}[0][2]][0] = $start;	# Start
			$bed_array{$chr}[$bed_array{$chr}[0][2]][1] = $end;	# End
			$bed_array{$chr}[0][2]++;
		}
	}
}
%bed_hash = (); # Empty the hash

close BED;
print "Finished loading input BED file.\n\n";

#########################################################################
#                         Reading Input Folders                         #
# For each chromosome found in the array, search for data files for the #
# chromosome in each folder. Load the chromosome input file from each   #
# folder and fill in %outhash with information from the bed file using  #
# the %bed_array: (1) chromosome (2) start (3) end & information from   #
# the input file: (4) header name (5) percentage methylation data.      #
# Print and clear %outhash for each chromosome, taking into account the #
# thresholds set by the user.                                           #
#########################################################################

# Print header
print OUT "Chromosome\tStart\tEnd";
for(my $i = 0; $i < @headernames; $i++){
	print OUT "\t$headernames[$i]";
}
print OUT "\n";

# Run process and print for each chromosome
foreach my $chr (sort keys %bed_array){
	print "Loading $chr\n";
	my %outhash;
	# For each sample/folder, run the chromosome bed file
	for(my $i = 0; $i < @filenames; $i++){
		my $filename = $filenames[$i] . $chr . ".bed";	# Create File Name
		my $arraycount = 0;	# Starting position for each iteration
		
		open (IN, "<$filename") or die "AvgMethv4.pl: Error: Couldn't open chromosome file $filename\n";
		my $inline = <IN>;  	# Check for header line
		if($inline =~ /^chr/){	# Process line if data
			my @inlinearray = split("\t",$inline);
			# Check the bed array for methylation information
			for(my $k = $arraycount; $k < $bed_array{$chr}[0][2]; $k++){
				# If end of in line is less than start of bed line, go to next in line
				if($inlinearray[2] < $bed_array{$chr}[$k][0]){last;}
				# If end of bed line is less than start of in line, go to next bed line
				# Increment starting position, previous lines no longer necessary for further in lines
				if($bed_array{$chr}[$k][1] < $inlinearray[1]){
					$arraycount++;
					next;
				}
				# If in line position lies within bed line position, add methylation information if sufficient
				if($inlinearray[1] >= $bed_array{$chr}[$k][0] && $inlinearray[2] <= $bed_array{$chr}[$k][1]){
					my @dash = split("-", $inlinearray[3]);
					# Check if above read threshold
					if ($dash[1] >= $minreads) {
						$outhash{$chr}{$bed_array{$chr}[$k][0]}{$bed_array{$chr}[$k][1]}{$filenames[$i]} .= ",$dash[0]";
					}
				}
				
			}
		}

		while(<IN>){
			my @inlinearray = split("\t",$_);
			# Check the bed array for methylation information
			for(my $k = $arraycount; $k < $bed_array{$chr}[0][2]; $k++){
				# If end of inline is less than start of bed line, go to next inline
				if($inlinearray[2] < $bed_array{$chr}[$k][0]){last;}
				# If end of bed line is less than start of in line, go to next bed line
				# Increment starting position, previous lines no longer necessary for further in lines
				if($bed_array{$chr}[$k][1] < $inlinearray[1]){
					$arraycount++;
					next;
				}
				# If in line position lies within bed line position, add methylation information if sufficient
				if($inlinearray[1] >= $bed_array{$chr}[$k][0] && $inlinearray[2] <= $bed_array{$chr}[$k][1]){
					my @dash = split("-", $inlinearray[3]);
					# Check if above read threshold
					if ($dash[1] >= $minreads) {
						$outhash{$chr}{$bed_array{$chr}[$k][0]}{$bed_array{$chr}[$k][1]}{$filenames[$i]} .= ",$dash[0]";
					}
				}
			}
		}
		close BED;

		# Check the entire data set for positions with no methylation information and fill in with "NA"
		for(my $k = 0; $k < $bed_array{$chr}[0][2]; $k++){
			if (!defined $outhash{$chr}{$bed_array{$chr}[$k][0]}{$bed_array{$chr}[$k][1]}{$filenames[$i]}){
				$outhash{$chr}{$bed_array{$chr}[$k][0]}{$bed_array{$chr}[$k][1]}{$filenames[$i]} = "NA";
			}
		}
	}

	print "Printing $chr\n";
	# Print data for this chromosome
	foreach my $outchr (keys %outhash){
		foreach my $outstart (sort {$a<=>$b} keys %{$outhash{$outchr}}){
			foreach my $outend (sort {$a<=>$b} keys %{$outhash{$outchr}{$outstart}}){
				# Make a count for file threshold
				my $maxNA = @filenames - $minfiles;
				for(my $i = 0; $i < @filenames; $i++) {
					# For each file missing data, increment count
					if($outhash{$outchr}{$outstart}{$outend}{$filenames[$i]} eq "NA") {
						$maxNA--;
					}
					# Files with data, process percentage methylation
					else {
						my @data = split(",",$outhash{$outchr}{$outstart}{$outend}{$filenames[$i]});
						# Check CpG site threshold
						if ($#data >= $mincpg){
							# Calculate percentage methylation
							my $sum = 0;
							foreach my $methyl (@data){
								if($methyl ne "") {$sum += $methyl;}
							}
							$outhash{$outchr}{$outstart}{$outend}{$filenames[$i]} = sprintf("%.5f", $sum/$#data);
						}
						# Does not pass threshold, add to missing data
						else {
							$outhash{$outchr}{$outstart}{$outend}{$filenames[$i]} = "NA";
							$maxNA--;
						}
					}
				}
				# Check file threshold, if met, then print
				if($maxNA >= 0) {
					# Print to out
					print OUT"$outchr\t$outstart\t$outend";
					for(my $i = 0; $i < @filenames; $i++) {
						print OUT "\t$outhash{$outchr}{$outstart}{$outend}{$filenames[$i]}";
					}
					print OUT "\n";
				}
			}
		}
	}


}

close OUT;
