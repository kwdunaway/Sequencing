#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu
# Date: 2-24-2014
# Script Name: AvgMeth.pl
#
# This script calculates average percent methylation of all CpG sites in each read of a 
# BED file.
#
# Arguments:
#    <see below>
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: AvgMeth.pl
    1) Output file
    2) BED or GTF file indicating areas of genome to average
    3) Minimum CpG Coung Threshold (must be at least 1)
    4+) Input Percent Methylation Folder (files need to have headers)
" unless @ARGV > 3;

my $outputname = shift(@ARGV);	# Output with average percentage methylation per PMD
open(OUT, ">$outputname") or die "AveragePercentMethylationAcrossBED.pl: Error: cannot open $outputname output file";
my $inputname = shift(@ARGV);	# Input PMD BED File
open(IN, "<$inputname") or die "AveragePercentMethylationAcrossBED.pl: Error: cannot open $inputname input BED file";
my $threshold = shift(@ARGV);	# Avg % Meth = "NA" unless >= minimum CpG site threshold
if($threshold < 1) {die "Threshold was less than 1";}

# Special case (default): 0 as threshold will cause 0 CpG sites to give "NA"
my @Folders = @ARGV;


# Global Variables
my @chrarray; 	   # Stores every field of lines in BED File and percentage methylation
my $linecount = 0; # Line count of the BED File

############################################
#           Reading Input BED File         #
############################################
# Inputs all of the data in the BED file 

my $firstline = <IN>;
if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
	print "No header found, processing first line of $inputname\n";
	my @line = split("\t",$firstline);
	$chrarray[$linecount][0][0] = $line[0];	# Chromosome
	$chrarray[$linecount][0][1] = $line[1];	# Start
	$chrarray[$linecount][0][2] = $line[2];	# End
	for(my $i = 0; $i < $#ARGV+1; $i++)		# For each input folder
	{
		# Initialize
		$chrarray[$linecount][$i+1][0] = "NA";	# Average percentage methylation
		$chrarray[$linecount][$i+1][1] = 0;	# Total Methylation
		$chrarray[$linecount][$i+1][2] = 0;	# CpG sites
	}
	$linecount++;
}
else { # If first line IS a header line
	print "Header Found in $inputname\n";
}

while(<IN>)
{ 
	chomp;
	my @line = split("\t",$_);
	
	$chrarray[$linecount][0][0] = $line[0];	# Chromosome
	$chrarray[$linecount][0][1] = $line[1];	# Start
	$chrarray[$linecount][0][2] = $line[2];	# End
	for(my $i = 0; $i < @Folders; $i++)		# For each input folder
	{
		# Initialize
		$chrarray[$linecount][$i+1][0] = "NA";	# Average percentage methylation
		$chrarray[$linecount][$i+1][1] = 0;	# Total Methylation
		$chrarray[$linecount][$i+1][2] = 0;	# CpG sites
	}
	$linecount++;
}
print "Loaded $linecount lines of $inputname \n";

############################################
#           Reading Input Folders          #
############################################


for(my $i = 0; $i < @Folders; $i++)	# Run process for each folder
{
	#           Obtain File Input Prefix

	my $inputprefix; # File prefix

	# Scan BED directory for number of chromosome files
#	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $Folders[$i];
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /(.+)chr./; 	# Extract file prefix
		$inputprefix = $1;
#		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
#		push (@Chr, $1); # Add to list of chromosome numbers
	}
#	@Chr = sort @Chr; # Sort list of chromosome numbers

	#                   Variables

	my $count = 0; # Count for BED reading loop
	my $filename; # Chromosome file name

	#           Process Input of Each File          

	while($count < $linecount)	# Run until end of BED file
	{
		# Open chromosome file
		
		#Make it work with capital chromosome names
		my $chrom = $chrarray[$count][0][0];
		$filename = $inputprefix . $chrom . ".bed";

		if (-e $filename) {
			open (BED, "<$filename") or die "Error: Couldn't open chromosome file $filename\n";
			print "Reading $filename \n";
#			print "Count: $count \t $chrarray[$count][0][0] \t $chrarray[$count][0][1] \t $chrarray[$count][0][2] \n";
		 } 
		else{
		#Gets rid of any sections of bed file that are outside the PerMeth file
			print "File does not exist: $filename \n";
			while($chrarray[$count][0][0] eq $chrom){
				$count++;
				if($count >= $linecount) {last;}
			}
			next;
		}


		my $bedline = <BED>;  	# Get rid of header line
		while(<BED>)
		{
			my @bedlinearray = split("\t",$_);
			# If undefined, end reading of this file
			if(!defined $bedlinearray[0] || !defined $chrarray[$count][0][0])
			{
				last;
			}
			# If diff chr, end reading of this file
			if($bedlinearray[0] ne $chrarray[$count][0][0])	
			{
				last;
			}
			# Go to next line if nothing found
			if($bedlinearray[2] > $chrarray[$count][0][2]) 
			{
				$chrarray[$count][$i+1][0] = "NA";
				$count++;
			} 
			else
			{
				# If within the range, then add to total meth
				if($bedlinearray[1] >= $chrarray[$count][0][1] && $bedlinearray[2] <= $chrarray[$count][0][2])
				{
					my @dash = split("-", $bedlinearray[3]);
					$chrarray[$count][$i+1][1] += $dash[0];	# Add percentage methylation
					$chrarray[$count][$i+1][2]++;			# Add to total sites
				}
				# Otherwise, skip to starting location of next region to get in the range
			}
		}

		#Gets rid of any sections of bed file that are outside the PerMeth file
		while($chrarray[$count][0][0] eq $chrom){
			$count++;
		}

		close BED;
	}
}

############################################
#     Calculate Percentage Methylation     #
############################################


# Prints header for results file
print OUT "chr\tstart\tend";
for(my $a = 0; $a < @Folders; $a++){
	print OUT "\t", $Folders[$a];
}
print OUT "\n";


#         Print results into OUT file
for(my $i = 0; $i < $linecount; $i++)
{
	print OUT $chrarray[$i][0][0], "\t", $chrarray[$i][0][1], "\t", $chrarray[$i][0][2];

	for(my $k = 1; $k <= @Folders; $k++)
	{
		# Sites under threshold -> "NA"
		if($chrarray[$i][$k][2] < $threshold)	
		{
			print OUT "\tNA";
		}
		# Sites Found -> Calculate average percent methylation
		else				
		{
			$chrarray[$i][$k][0] = ($chrarray[$i][$k][1])/($chrarray[$i][$k][2]); # Calculate average
			my $rounded = sprintf("%.5f", $chrarray[$i][$k][0]);
			print OUT "\t", $rounded;
		}
	}
	print OUT "\n";
}

close IN;
close OUT;


__END__

