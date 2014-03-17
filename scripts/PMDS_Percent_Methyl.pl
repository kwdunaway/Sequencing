#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Roy Chu
# Email: rgchu@ucdavis.edu
# Date: 12-3-2013
# Script Name: PMDS_Percent_Methyl.pl
#
# This script takes a bed file and looks at % methylation across each region in the bed file.
#
# Arguments:
#    1) Input BED file
#    2) Input Percent_Methyl Folder
#    3) Output file
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: PMDS_Percent_Methyl.pl 
    1) Input BED file with PMDs
    2) Input Percent_Methyl Folder
    3) Output file
" unless @ARGV == 3;

my $inputname = shift(@ARGV);	# Input PMD BED File
open(IN, "<$inputname") or die "PMDS_Percent_Methyl.pl: Error: cannot open $inputname GTF infile";
my $inputfolder = shift(@ARGV); # Percent Methyl Folder
my $outputname = shift(@ARGV);	# Output with average percentage methylation per PMD
open(OUT, ">$outputname") or die "PMDS_Percent_Methyl.pl: Error: cannot open $outputname OUT outfile";

my $inputprefix; # File prefix
# Scan BED directory for number of chromosome files
my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
my $filedir = $inputfolder;
$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
my @files = glob( $filedir . '*' ); # Gets list of all files in directory
@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

foreach my $file (@files) {
	$file =~ /(.+)chr./; 	# Extract file prefix
	$inputprefix = $1;
	$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
	push (@Chr, $1); # Add to list of chromosome numbers
}
@Chr = sort @Chr; # Sort list of chromosome numbers

# Global Variables
my @chrarray; # Stores all PMD information
my $linecount = 0; # Count for each location
my $count = 0; # Count for BED reading loop
my $filename; # Chromosome file name



##############################
# 	      Header             #
##############################

my $header = <IN>;
if ($header =~ /chr/)	# No header, first line is data
{
	print OUT "Chromosome\tStart\tEnd\n";
}
else			# Header found
{
	print OUT $header;
	$header = <IN>;
	my @line = split("\t",$header);
	$chrarray[$linecount][0] = $line[0];	# Chromosome
	$chrarray[$linecount][1] = $line[1];	# Start
	$chrarray[$linecount][2] = $line[2];	# End	
	$chrarray[$linecount][3] = "NA";	# Average percentage methylation
	$chrarray[$linecount][4] = 0;		# Total Methylation
	$chrarray[$linecount][5] = 0;		# CpG sites
	$linecount++;
}



##############################
# 	   Store Input           #
##############################

# Get first line of data from IN
while(<IN>){ # Store input file in 2D array
	chomp;
	my @line = split("\t",$_);
	$chrarray[$linecount][0] = $line[0];	# Chromosome
	$chrarray[$linecount][1] = $line[1];	# Start
	$chrarray[$linecount][2] = $line[2];	# End
	$chrarray[$linecount][3] = "NA";	# Average percentage methylation
	$chrarray[$linecount][4] = 0;		# Total Methylation
	$chrarray[$linecount][5] = 0;		# CpG sites
	$linecount++;
}



##############################
# 	  Process Input          #
##############################

# Read all chromosomes found in PMD file
while($count < $linecount)
{
	$filename = $inputprefix . $chrarray[$count][0] . ".bed";

	open (BED, "<$filename") or die "PMDS_Percent_Methyl.pl: Error: Couldn't open chromosome file $filename\n";

	# Scan folder
	my $bedline = <BED>;  	# Get rid of header line
	while(<BED>)
	{
		my @bedlinearray = split("\t",$_);
		if(!defined $bedlinearray[0] || !defined $chrarray[$count][0])
		{
			last;
		}
		if($bedlinearray[0] ne $chrarray[$count][0])	# If diff chr, end reading
		{
			last;
		}
		if($bedlinearray[2] > $chrarray[$count][2]) # Go to next PMD if nothing found
		{
			$chrarray[$count][3] = "NA";
			$count++;
		} 
		else
		{
			# If within the range, then add to total meth
			if($bedlinearray[1] >= $chrarray[$count][1] && $bedlinearray[2] <= $chrarray[$count][2])
			{
				my @dash = split("-", $bedlinearray[3]);
				$chrarray[$count][4] += $dash[0];	# Add percentage methylation
				$chrarray[$count][5]++;			# Add to total sites
			}
			# Otherwise, skip to starting location of PMD to get in the range
		}
	}
	close BED;
}

# Print results into OUT file
for(my $i = 0; $i < $count; $i++)
{
	if($chrarray[$i][5] == 0)	# No Sites, "NA"
	{
		print OUT $chrarray[$i][0], "\t", $chrarray[$i][1], "\t", $chrarray[$i][2], "\t", $chrarray[$i][3], "\n";
	}
	else				# Sites Found
	{
		$chrarray[$i][3] = ($chrarray[$i][4])/($chrarray[$i][5]); # Calculate average
		print OUT $chrarray[$i][0], "\t", $chrarray[$i][1], "\t", $chrarray[$i][2], "\t", $chrarray[$i][3], "\n";
	}
}

close IN;
close OUT;


__END__

