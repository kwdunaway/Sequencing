#!/usr/bin/perl 
use strict; use warnings;

################################################################################
# Author: Keith Dunaway and Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 8-14-2014
#
# This script calculates average percent methylation of all CpG sites in each 
# line of a BED file. Multiple percent methylation folders with 
# percent methylationbed files of each chromosome may be entered as inputs to
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
my %bed_hash;	# Stores every field of lines in BED File and percentage methylation

############################################
#           Reading Input BED File         #
############################################
# Inputs all of the data in the BED file 

my $firstline = <IN>;		# Check for header
if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
	print "No header found, processing first line.\n";
	chomp($firstline);
	my @line = split("\t",$firstline);
	if(!defined $bed_hash{$line[0]}{$line[1]}{$line[2]}) {
		$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Chromosome
	}
	else {
		print "Duplicate found, skipping: $firstline";
	}
}
else { # If first line IS a header line
	print "Header found, skipping first line!\n";
}

# Process entire BED file
while(<IN>)
{ 
	chomp;
	my @line = split("\t",$_);
	if ($line[0] =~ /_/){next;}
	if(!defined $bed_hash{$line[0]}{$line[1]}{$line[2]}) {
		$bed_hash{$line[0]}{$line[1]}{$line[2]} = 1;	# Chromosome
	}
	else {
		print "Duplicate found, skipping: $_";
	}
}
close IN;
print "past loading input BED file\n";

############################################
#           Reading Input Folders          #
############################################

# Print header
print OUT "Chromosome\tStart\tEnd";
for(my $i = 0; $i < $#ARGV+1; $i++){
	print OUT "\t$ARGV[$i]";
}
print OUT "\n";

# Run process and print for each chromosome
foreach my $chr (sort keys %bed_hash){
	print "Loading $chr\n";
	my %outhash;
	foreach my $start (keys %{$bed_hash{$chr}}){
		foreach my $end (keys %{$bed_hash{$chr}{$start}}){
			my $entry = $chr . ":" . $start . "-". $end;
			for(my $i = 0; $i < $#ARGV+1; $i++){	# Run line for each folder
				# Open chromosome file
				my $filename = $ARGV[$i] . $chr . ".bed";

				open (BED, "<$filename") or die "AvgMethv4.pl: Error: Couldn't open chromosome file $filename\n";
		
				my $bedline = <BED>;  	# Check for header line
				if($bedline =~ /^chr/){	# Process line if data
					my @bedlinearray = split("\t",$bedline);
					# Go to next line if nothing found
					if($bedlinearray[1] >= $start && $bedlinearray[2] <= $end){
						my @dash = split("-", $bedlinearray[3]);
						if ($dash[1] >= $minreads) {
							$outhash{$entry}{$ARGV[$i]} .= ",$dash[0]";
						}
					}
				}

				while(<BED>){
					my @bedlinearray = split("\t",$_);
					# Go to next line if nothing found
					if($bedlinearray[1] >= $start && $bedlinearray[2] <= $end){
						my @dash = split("-", $bedlinearray[3]);
						if ($dash[1] >= $minreads) {
							$outhash{$entry}{$ARGV[$i]} .= ",$dash[0]";
						}
					}
				}
				close BED;

				if (!defined $outhash{$entry}{$ARGV[$i]}){
					$outhash{$entry}{$ARGV[$i]} = "NA";
				}
			}
		}
	}

	print "Printing $chr\n";
	# Print data for this chromosome
	foreach my $entry (keys %outhash){
		my $maxNA = $#ARGV + 1 - $minfiles;
		for(my $i = 0; $i < $#ARGV+1; $i++) {
			if($outhash{$entry}{$ARGV[$i]} eq "NA") {
				$maxNA--;
			}
			else {
				my @data = split(",",$outhash{$entry}{$ARGV[$i]});
				if ($#data >= $mincpg){#Check CpG site threshold
					my $sum = 0;
					foreach my $methyl (@data){
						if($methyl ne "") {$sum += $methyl;}
					}
					$outhash{$entry}{$ARGV[$i]} = sprintf("%.5f", $sum/$#data);
				}
				else {
					$outhash{$entry}{$ARGV[$i]} = "NA";
					$maxNA--;
				}
			}
		}
		if($maxNA >= 0) {
			# Print to out
			my @chr = split(":",$entry);
			my @pos = split("-",$chr[1]);
			print OUT"$chr[0]\t$pos[0]\t$pos[1]";
			for(my $i = 0; $i < $#ARGV+1; $i++) {
				print OUT "\t$outhash{$entry}{$ARGV[$i]}";
			}
			print OUT "\n";
		}
	}


}

close OUT;
