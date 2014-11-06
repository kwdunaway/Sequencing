#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu kwdunaway@ucdavis.edu
# Date: 11-5-2014
#
# This script counts CpG sites in areas of the genome specified by an input 
# file (with columns chr, start, end). An input file with the genome will be
# scanned for CpG sites in these regions and CpG islands can be masked out.
# If the masking option is chosen, there will be a fifth column describing
# the number of bases in the region that are unmasked.
#
# Arguments:
#    <see below>
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: cpgcounter.pl
    1) Input File
    2) Genome fasta file (Example: hg19.fa)
    3) CpG Islands mask file (Optional, type 'NA' if none)
    4) Output File
" unless @ARGV == 4;

my $inputfile = shift(@ARGV);	# Input file with tab-delimited columns(chromosome, start, end)
my $fasta = shift(@ARGV);	# FASTA file containing a genome
my $cpgislands = shift(@ARGV);	# File containing CpG islands to be masked (same format as input)
my $outputfile = shift(@ARGV);	# Output file name
my $cg = "CG";


open(IN, "<$inputfile") or die "cpgcounter.pl: Error: cannot open $inputfile input file";
open(GEN, "<$fasta") or die "cpgcounter.pl: Error: cannot open $fasta genome file";
open(OUT, ">$outputfile") or die "cpgcounter.pl: Error: cannot open $outputfile output file";

print "Initializing\n";

my %areahash; # Hash to store the input file information

# Procedure for No CpG Islands Mask File
if($cpgislands eq "NA") {
	my $firstline = <IN>;		# Check for header
	if ($firstline =~ /^chr[0-9]/){	# Checks to see if the first line is not a header
		print "No header found in $inputfile, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {} # Ignore non-standard chromosomes
		else { # Push line to hash
			$areahash{$line[0]}[0][2] = 0;	# Initialize count of this chromosome
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1]; # Start
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2]; # End
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0; # CpG site count
			$areahash{$line[0]}[0][2]++; # Increment count for next line
		}	
	}
	else { # If first line IS a header line
		print "Header found in $inputfile, skipping first line.\n";
	}

	print "Reading input areas.\n";
	while(<IN>){
		chomp;
		my @line = split("\t",$_);
		if ($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
		# Initialize count of this chromosome if not defined already
		if(!defined $areahash{$line[0]}[0][2]) {$areahash{$line[0]}[0][2] = 0;}
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1]; # Start
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2]; # End
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0; # CpG site count
		$areahash{$line[0]}[0][2]++; # Increment count for next line
	}
	close IN;

	print "Scanning Genome for CpG Sites.\n";
	my $currchr;		# Iterator for current chromosome
	my $count = 0;		# Region count for input hash
	my $lowerend = 0;	# Lower coordinate of line
	my $higherend = 0;	# Upper coordinate of line
	my $prevchar = "Z";	# Last character of previous line

	# Reading genome file
	while(<GEN>){
		chomp;
		my $line = $_;
		if($line =~ /_/) {next;} # Ignore non-standard chromosomes
		if($line =~ />chr/) { # New chromosome
			# Re-initialize variables
			$line =~ />(chr.+)/;
			$currchr = $1;
			$count = 0;
			$lowerend = 0;
			$higherend = 0;
			$prevchar = "Z";
			print "Scanning ", $currchr, "\n";
			next;
		}
		# Chromosome in input not found in genome
		if(!defined $areahash{$currchr}[0][2]) {next;} 
		$line = uc($line); # Uppercase line for consistency
		# Case 1: Line is after region, move to next region
		while($count < $areahash{$currchr}[0][2] && $lowerend > $areahash{$currchr}[$count][1]){
			$count++;
		}
		# No more input regions for this chromosome
		if($count == $areahash{$currchr}[0][2]) {next;}
		# Increment higher end
		$higherend = $higherend + length($line);
		# Case 2: Line is before region, move to next line
		if($higherend < $areahash{$currchr}[$count][0]) {
			# Increment lower end
			$lowerend = $lowerend + length($line);
			# Set new previous character
			$prevchar = substr($line, -1);
			next;
		}
		# Case 3: Line encompasses region
		if($lowerend <= $areahash{$currchr}[$count][0] && $higherend >= $areahash{$currchr}[$count][1]){
			my $low = $areahash{$currchr}[$count][0] - $lowerend;
			my $high = $areahash{$currchr}[$count][1] - $areahash{$currchr}[$count][0];
			my $area = substr($line, $low, $high); # Take enclosed area
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
		}
		# Case 4: Region encompasses line
		elsif($lowerend >= $areahash{$currchr}[$count][0] && $higherend <= $areahash{$currchr}[$count][1]){
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $line =~ /$cg/g;
			# Check if edge may contain CG
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
		}
		# Case 5: The higher end of the line overlaps with the lower end of the region
		elsif($higherend > $areahash{$currchr}[$count][0] && $higherend < $areahash{$currchr}[$count][1]){
			my $high = -($higherend - $areahash{$currchr}[$count][0]);
			my $area = substr($line, $high);
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;	
		}
		# Case 6: The higher end of the region overlaps with the lower end of the line
		elsif($lowerend > $areahash{$currchr}[$count][0] && $lowerend < $areahash{$currchr}[$count][1]){
			my $low = $areahash{$currchr}[$count][1] - $lowerend;
			my $area = substr($line, 0, $low);
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			# Check if edge may contain CG
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
		}
		# Increment lower end
		$lowerend = $lowerend + length($line);
		# Set new previous character
		$prevchar = substr($line, -1);
	}
	close GEN;

	print "Printing to output.\n";
	foreach my $allchr (sort keys %areahash){ # Sort each chromosome
		# Skip undefined chromosomes
		if(!defined $areahash{$allchr}[0][2]) {next;}
		for(my $k = 0; $k < $areahash{$allchr}[0][2]; $k++){
			print OUT $allchr, "\t", $areahash{$allchr}[$k][0], "\t", $areahash{$allchr}[$k][1], "\t", $areahash{$allchr}[$k][3], "\n";
		}
	}
}
# Procedure for Masking CpG Islands
else {
	my %cgi_unsorted; # Unsorted hash to take in information from CpG islands
	my %cgi_sorted; # Sorted version of the CpG island hash
	open(CGI, "<$cpgislands") or die "cpgcounter.pl: Error: cannot open $cpgislands CpG Islands file";
	# Check for header
	my $firstline = <CGI>;		
	if ($firstline =~ /^chr[0-9]/){	# Checks to see if the first line is not a header
		print "No header found in $cpgislands, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {} # Ignore non-standard chromosomes
		else {$cgi_unsorted{$line[0]}{$line[1]}{$line[2]} = 1;}	# Push line to hash
	}
	else { # If first line IS a header line
		print "Header found in $cpgislands, skipping first line.\n";
	}

	# Process Masking File
	print "Reading CpG islands file.\n";
	while(<CGI>){
		chomp;
		my @line = split("\t",$_);
		if($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
		# Check if duplicate, otherwise enter into hash
		if(!defined $cgi_unsorted{$line[0]}{$line[1]}{$line[2]}) {
			$cgi_unsorted{$line[0]}{$line[1]}{$line[2]} = 1;
		}
		# Otherwise duplicate
		else {print "Warning: duplicate found, skipping: $_\n";}
	}
	close CGI;

	# Fill in hash-array with hash information, this is done to sort the information
	foreach my $chr (keys %cgi_unsorted){
		# Initialize line count, varies for each chromosome
		$cgi_sorted{$chr}[0][2] = 0;
		# For each position, sort numerically
		foreach my $start (sort {$a<=>$b} keys %{$cgi_unsorted{$chr}}){
			foreach my $end (sort {$a<=>$b} keys %{$cgi_unsorted{$chr}{$start}}){
				$cgi_sorted{$chr}[$cgi_sorted{$chr}[0][2]][0] = $start;	# Start
				$cgi_sorted{$chr}[$cgi_sorted{$chr}[0][2]][1] = $end;	# End
				$cgi_sorted{$chr}[0][2]++; # Increment count
			}
		}
	}
	%cgi_unsorted = (); # Empty the unsorted hash

	# Reading input file
	$firstline = <IN>;		# Check for header
	if ($firstline =~ /^chr[0-9]/){	# Checks to see if the first line is not a header
		print "No header found in $inputfile, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {} # Ignore non-standard chromosomes
		else { # Push line to hash
			$areahash{$line[0]}[0][2] = 0; # Initialize line count of this chromosome
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1]; # Start
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2]; # End
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0; # CpG site count 
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][4]=0; # Count for unmasked bases
			$areahash{$line[0]}[0][2]++; # Increment line(region) count
		}	
	}
	else { # If first line IS a header line
		print "Header found in $inputfile, skipping first line.\n";
	}

	print "Reading input areas file.\n";
	while(<IN>){
		chomp;
		my @line = split("\t",$_);
		if ($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
		# Initialize count of this chromosome if not defined already
		if(!defined $areahash{$line[0]}[0][2]) {$areahash{$line[0]}[0][2] = 0;}
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1]; # Start
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2]; # End
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0; # CpG site count 
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][4]=0; # Count for unmasked bases
		$areahash{$line[0]}[0][2]++; # Increment line(region) count
	}
	close IN;

	print "Scanning Genome for CpG Sites.\n";
	my $currchr;		# Iterator for current chromosome
	my $count = 0;		# Region count for input hash
	my $lowerend = 0;	# Lower coordinate of line
	my $higherend = 0;	# Upper coordinate of line
	my $prevchar = "Z";	# Last character of previous line
	my $cpgpositioncounter = 0;	# Increments starting position for cpg island checking

	# Reading genome file
	while(<GEN>){
		chomp;
		my $line = $_;
		if($line =~ /_/) {next;} # Ignore non-standard chromosomes
		if($line =~ />chr/) {	# New chromosome
			# Re-initialize variables
			$line =~ />(chr.+)/;
			$currchr = $1;
			$count = 0;
			$lowerend = 0;
			$higherend = 0;
			$prevchar = "Z";
			$cpgpositioncounter = 0;
			print "Scanning ", $currchr, "\n";
			next;
		}
		# Chromosome in input not found in genome
		if(!defined $areahash{$currchr}[0][2]) {next;}
		$line = uc($line); # Uppercase line for consistency
		# Case 1: Line is after region, move to next region
		while($count < $areahash{$currchr}[0][2] && $lowerend > $areahash{$currchr}[$count][1]){
			$count++;
		}
		# No more input regions for this chromosome
		if($count == $areahash{$currchr}[0][2]) {next;}
		# Increment higher end
		$higherend = $higherend + length($line);
		# Case 2: Line is before region, move to next line
		if($higherend < $areahash{$currchr}[$count][0]) {
			# Increment lower end
			$lowerend = $lowerend + length($line);
			# Set new previous character
			$prevchar = substr($line, -1);
			next;
		}
		# Case 3: Line encompasses region
		if($lowerend <= $areahash{$currchr}[$count][0] && $higherend >= $areahash{$currchr}[$count][1]){
			# Check for CpG Islands, Transliterate CGs into Ns to mask them
			for(my $k = $cpgpositioncounter; $k < $cgi_sorted{$currchr}[0][2]; $k++){
				# Case 1: CpG Island is after region, end CpG loop
				if($cgi_sorted{$currchr}[$k][0] > $higherend){last;}
				# Case 2: CpG Island is before region, go to next CpG Island
				if($cgi_sorted{$currchr}[$k][1] < $lowerend){
					$cpgpositioncounter++;
					next;
				}
				# Case 3: Line encompasses CpG Island
				if($lowerend <= $cgi_sorted{$currchr}[$k][0] && $higherend >= $cgi_sorted{$currchr}[$k][1]){
					my $low = $cgi_sorted{$currchr}[$k][0] - $lowerend;
					my $high = $cgi_sorted{$currchr}[$k][1] - $cgi_sorted{$currchr}[$k][0];
					substr($line,$low,$high) =~ tr/ATCG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/ATCG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/ATCG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/ATCG/N/;
				}
			}
			my $low = $areahash{$currchr}[$count][0] - $lowerend;
			my $high = $areahash{$currchr}[$count][1] - $areahash{$currchr}[$count][0];
			my $area = substr($line, $low, $high);
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			# Count number of masked bases
			$areahash{$currchr}[$count][4] += () = $area =~ /N/g;
		}
		# Case 4: Region encompasses line
		elsif($lowerend >= $areahash{$currchr}[$count][0] && $higherend <= $areahash{$currchr}[$count][1]){
			# Check for CpG Islands, Transliterate CGs into Ns to mask them
			for(my $k = $cpgpositioncounter; $k < $cgi_sorted{$currchr}[0][2]; $k++){
				# Case 1: CpG Island is after region, end CpG loop
				if($cgi_sorted{$currchr}[$k][0] > $higherend){last;}
				# Case 2: CpG Island is before region, go to next CpG Island
				if($cgi_sorted{$currchr}[$k][1] < $lowerend){
					$cpgpositioncounter++;
					next;
				}
				# Case 3: Line encompasses CpG Island
				if($lowerend <= $cgi_sorted{$currchr}[$k][0] && $higherend >= $cgi_sorted{$currchr}[$k][1]){
					my $low = $cgi_sorted{$currchr}[$k][0] - $lowerend;
					my $high = $cgi_sorted{$currchr}[$k][1] - $cgi_sorted{$currchr}[$k][0];
					substr($line,$low,$high) =~ tr/ATCG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/ATCG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/ATCG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/ATCG/N/;
				}
			}
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $line =~ /$cg/g;
			# Count number of masked bases
			$areahash{$currchr}[$count][4] += () = $line =~ /N/g;
			# Check edge case for CG
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
		}
		# Case 5: The higher end of the line overlaps with the lower end of the region
		elsif($higherend > $areahash{$currchr}[$count][0] && $higherend < $areahash{$currchr}[$count][1]){
			# Check for CpG Islands, Transliterate CGs into Ns to mask them
			for(my $k = $cpgpositioncounter; $k < $cgi_sorted{$currchr}[0][2]; $k++){
				# Case 1: CpG Island is after region, end CpG loop
				if($cgi_sorted{$currchr}[$k][0] > $higherend){last;}
				# Case 2: CpG Island is before region, go to next CpG Island
				if($cgi_sorted{$currchr}[$k][1] < $lowerend){
					$cpgpositioncounter++;
					next;
				}
				# Case 3: Line encompasses CpG Island
				if($lowerend <= $cgi_sorted{$currchr}[$k][0] && $higherend >= $cgi_sorted{$currchr}[$k][1]){
					my $low = $cgi_sorted{$currchr}[$k][0] - $lowerend;
					my $high = $cgi_sorted{$currchr}[$k][1] - $cgi_sorted{$currchr}[$k][0];
					substr($line,$low,$high) =~ tr/ATCG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/ATCG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/ATCG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/ATCG/N/;
				}
			}
			my $high = -($higherend - $areahash{$currchr}[$count][0]);
			my $area = substr($line, $high);
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			# Count number of masked bases
			$areahash{$currchr}[$count][4] += () = $area =~ /N/g;		
		}
		# Case 6: The higher end of the region overlaps with the lower end of the line
		elsif($lowerend > $areahash{$currchr}[$count][0] && $lowerend < $areahash{$currchr}[$count][1]){
			# Check for CpG Islands, Transliterate CGs into Ns to mask them
			for(my $k = $cpgpositioncounter; $k < $cgi_sorted{$currchr}[0][2]; $k++){
				# Case 1: CpG Island is after region, end CpG loop
				if($cgi_sorted{$currchr}[$k][0] > $higherend){last;}
				# Case 2: CpG Island is before region, go to next CpG Island
				if($cgi_sorted{$currchr}[$k][1] < $lowerend){
					$cpgpositioncounter++;
					next;
				}
				# Case 3: Line encompasses CpG Island
				if($lowerend <= $cgi_sorted{$currchr}[$k][0] && $higherend >= $cgi_sorted{$currchr}[$k][1]){
					my $low = $cgi_sorted{$currchr}[$k][0] - $lowerend;
					my $high = $cgi_sorted{$currchr}[$k][1] - $cgi_sorted{$currchr}[$k][0];
					substr($line,$low,$high) =~ tr/ATCG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/ATCG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/ATCG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/ATCG/N/;
				}
			}
			my $low = $areahash{$currchr}[$count][1] - $lowerend;
			my $area = substr($line, 0, $low);
			# Count CpG Sites
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			# Count number of masked bases
			$areahash{$currchr}[$count][4] += () = $area =~ /N/g;
			# Check edge case for CG
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
		}
		# Increment lower end
		$lowerend = $lowerend + length($line);
		# Set new previous character
		$prevchar = substr($line, -1);
	}
	close GEN;

	print "Printing to output.\n";
	foreach my $allchr (sort keys %areahash){ # Sort chromosomes
		# Skip undefined chromosomes	
		if(!defined $areahash{$allchr}[0][2]) {next;} 
		for(my $k = 0; $k < $areahash{$allchr}[0][2]; $k++){
			# Calculate number of unmasked bases
			$areahash{$allchr}[$k][4] = ($areahash{$allchr}[$k][1] - $areahash{$allchr}[$k][0]) - $areahash{$allchr}[$k][4];
			print OUT $allchr, "\t", $areahash{$allchr}[$k][0], "\t", $areahash{$allchr}[$k][1], "\t", $areahash{$allchr}[$k][3], "\t", $areahash{$allchr}[$k][4], "\n";
		}
	}
}

close OUT;
