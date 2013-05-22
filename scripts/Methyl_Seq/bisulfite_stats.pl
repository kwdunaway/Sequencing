#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 7-31-2012
# Script Name: bisulfite_stats.pl
#
#
#
# Arguments:
#    See Below
#
################################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Input file name
    2) Output prefix
" unless @ARGV == 2;

my $inputfile = shift(@ARGV);
my $outputprefix = shift(@ARGV);

####################
#    Processing    #
####################

bisulfite_stats($inputfile,$outputprefix);

sub bisulfite_stats
{
	# Input
	my ($inputfile,$outputprefix) = @_;

	# File I/O
	my $temp_CG = $outputprefix . "_CGtemp.txt";
	my $temp_CHG = $outputprefix . "_CHGtemp.txt";
	my $temp_CHH = $outputprefix . "_CHHtemp.txt";
	my $CGout = $outputprefix . "_CG.txt";
	my $CHGout = $outputprefix . "_CHG.txt";
	my $CHHout = $outputprefix . "_CHH.txt";


	open(IN, "<$inputfile") or die "Error: Bisulfite Stats: cannot open $inputfile input file";
	open(CGOUT, ">$temp_CG") or die "Error: Bisulfite Stats: cannot open $temp_CG output file";
	open(CHGOUT, ">$temp_CHG") or die "Error: Bisulfite Stats: cannot open $temp_CHG output file";
	open(CHHOUT, ">$temp_CHH") or die "Error: Bisulfite Stats: cannot open $temp_CHH output file";

	# Hash{CG/CHG/CHH}->{Number of Sites}->[Defined Flag, Number of Reads, Number of Sites]
	# Each element in the array: 
		# [0] Flag to check if this number of sites has a defined array
		# [1] Number of reads with this number of sites
		# [2-# of Sites] Count of reads with [x-2] number of methylated sites, [2] is 0 methylated, [3] is 1 methylated, etc.

	my %CGhash; # Hash to hold CpG sites data
	my %CHGhash; # Hash to hold CHG sites data
	my %CHHhash; # Hash to hold CHH sites data

	while(<IN>){
		chomp;
		my @line = split("\t",$_);
		my $sites_line = $line[6];
		my $CGsites = 0; # Number of CG sites
		my $CHGsites = 0; # Number of CHG sites
		my $CHHsites = 0; # Number of CHH sites
		my $num_methylated = 2; # Counts number of methylated sites (Start at array number 2)
	
		# CG
		$CGsites++ while ($sites_line =~ m/[xX]/g); # Match sites
		$num_methylated++ while ($sites_line =~ m/X/g); # Match methylated sites
		if(! defined $CGhash{$CGsites}[0]) { # Check if array has been defined yet
			for(my $i = 0; $i <= $CGsites; $i++) {
				$CGhash{$CGsites}[$i+2] = 0;
			}
			$CGhash{$CGsites}[0] = 1; # Flag that array is now defined
		}
		$CGhash{$CGsites}[1]++; # Add to total reads
		$CGhash{$CGsites}[$num_methylated]++; # Add to read with that # of methylated sites

		# CHG
		$num_methylated = 2; # Reset to 2
		$CHGsites++ while ($sites_line =~ m/[yY]/g); # Match sites
		$num_methylated++ while ($sites_line =~ m/Y/g); # Match methylated sites
		if(! defined $CHGhash{$CHGsites}[0]) { # Check if array has been defined yet
			for(my $i = 0; $i <= $CHGsites; $i++) {
				$CHGhash{$CHGsites}[$i+2] = 0;
			}
			$CHGhash{$CHGsites}[0] = 1; # Flag that array is now defined
		}
		$CHGhash{$CHGsites}[1]++; # Add to total reads
		$CHGhash{$CHGsites}[$num_methylated]++; # Add to read with that # of methylated sites

		# CHH
		$num_methylated = 2; # Reset to 2
		$CHHsites++ while ($sites_line =~ m/[zZ]/g); # Match sites
		$num_methylated++ while ($sites_line =~ m/Z/g); # Match methylated sites
		if(! defined $CHHhash{$CHHsites}[0]) { # Check if array has been defined yet
			for(my $i = 0; $i <= $CHHsites; $i++) {
				$CHHhash{$CHHsites}[$i+2] = 0;
			}
			$CHHhash{$CHHsites}[0] = 1; # Flag that array is now defined
		}
		$CHHhash{$CHHsites}[1]++; # Add to total reads
		$CHHhash{$CHHsites}[$num_methylated]++; # Add to read with that # of methylated sites
	}

	# Description of each tab-delimited column
	# # of Sites	Total Reads	Count of Methylated Sites

	foreach my $sites (keys %CGhash) {
		print CGOUT $sites, "\t", $CGhash{$sites}[1], "\t";
		for (my $i = 0; $i <= $sites; $i++) {
			if($i != 0) {print CGOUT ",";}
			print CGOUT $CGhash{$sites}[$i+2];
		}
		print CGOUT "\n";
	}

	foreach my $sites (keys %CHGhash) {
		print CHGOUT $sites, "\t", $CHGhash{$sites}[1], "\t";
		for (my $i = 0; $i <= $sites; $i++) {
			if($i != 0) {print CHGOUT ",";}
			print CHGOUT $CHGhash{$sites}[$i+2];
		}
		print CHGOUT "\n";
	}

	foreach my $sites (keys %CHHhash) {
		print CHHOUT $sites, "\t", $CHHhash{$sites}[1], "\t";
		for (my $i = 0; $i <= $sites; $i++) {
			if($i != 0) {print CHHOUT ",";}
			print CHHOUT $CHHhash{$sites}[$i+2];
		}
		print CHHOUT "\n";
	}

	close IN;
	close CGOUT;
	close CHGOUT;
	close CHHOUT;

	# Sorting files and removing the temporary files

	`sort -nk 1 $temp_CG > $CGout`;
	`sort -nk 1 $temp_CHG > $CHGout`;
	`sort -nk 1 $temp_CHH > $CHHout`;
	`rm $temp_CG $temp_CHG $temp_CHH`;
}

__END__
