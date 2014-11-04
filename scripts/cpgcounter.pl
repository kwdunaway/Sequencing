#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu
# Date: 2-24-2014
# Script Name: genometofastq.pl
#
# This script counts CpG sites in areas of the genome specified by an input 
# file (with columns chr, start, end). An input file with the genome will be
# scanned for CpG sites in these regions and CpG islands can be masked out.
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

my $inputfile = shift(@ARGV);
my $fasta = shift(@ARGV);
my $cpgislands = shift(@ARGV);	
my $outputfile = shift(@ARGV);
my $cg = "CG";


open(IN, "<$inputfile") or die "cpgcounter.pl: Error: cannot open $inputfile input file";
open(GEN, "<$fasta") or die "cpgcounter.pl: Error: cannot open $fasta genome file";
open(OUT, ">$outputfile") or die "cpgcounter.pl: Error: cannot open $outputfile output file";

print "Initializing\n";

my %areahash;

# Procedure for No CpG Islands Mask File
if($cpgislands eq "NA") {
	my $firstline = <IN>;		# Check for header
	if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
		print "No header found in $inputfile, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {} # Ignore non-standard chromosomes
		else { # Push line to hash
			$areahash{$line[0]}[0][2] = 0;
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1];
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2];
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0;
			#print OUT $cgi{$line[0]}[$cgi{$line[0]}[0][2]][0], "\t", $cgi{$line[0]}[$cgi{$line[0]}[0][2]][1], "\n";
			$areahash{$line[0]}[0][2]++;
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
		if(!defined $areahash{$line[0]}[0][2]) {$areahash{$line[0]}[0][2] = 0;}
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1];
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2];
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0;
		#print OUT $cgi{$line[0]}[$cgi{$line[0]}[0][2]][0], "\t", $cgi{$line[0]}[$cgi{$line[0]}[0][2]][1], "\n";
		$areahash{$line[0]}[0][2]++;
	}
	close IN;

	print "Scanning Genome for CpG Sites.\n";
	my $currchr;
	my $count = 0;
	my $lowerend = 0;
	my $higherend = 0;
	my $prevchar = "Z";

	while(<GEN>){
		chomp;
		my $line = $_;
		if($line =~ /_/) {next;}
		if($line =~ />chr/) {	# New chromosome
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
		if(!defined $areahash{$currchr}[0][2]) {next;}
		$line = uc($line);
		while($count < $areahash{$currchr}[0][2] && $lowerend > $areahash{$currchr}[$count][1]){
			$count++;
		}
		if($count == $areahash{$currchr}[0][2]) {next;}
		$higherend = $higherend + length($line);
		if($higherend < $areahash{$currchr}[$count][0]) {
			$lowerend = $lowerend + length($line);
			$prevchar = substr($line, -1);
			next;
		}
		if($lowerend <= $areahash{$currchr}[$count][0] && $higherend >= $areahash{$currchr}[$count][1]){
			my $low = $areahash{$currchr}[$count][0] - $lowerend;
			my $high = $areahash{$currchr}[$count][1] - $areahash{$currchr}[$count][0];
			my $area = substr($line, $low, $high);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";
		}
		elsif($lowerend >= $areahash{$currchr}[$count][0] && $higherend <= $areahash{$currchr}[$count][1]){
			#print "Squish\n";
			$areahash{$currchr}[$count][3] += () = $line =~ /$cg/g;
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\n";	
			#print $line, "\t", $areahash{$currchr}[$count][3], "\n";
		}
		elsif($higherend > $areahash{$currchr}[$count][0] && $higherend < $areahash{$currchr}[$count][1]){
			#print "High\n";
			my $high = -($higherend - $areahash{$currchr}[$count][0]);
			my $area = substr($line, $high);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {$areahash{$currchr}[$count][3]++;}
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\t", "\t", $high, "\n";
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";		
		}
		elsif($lowerend > $areahash{$currchr}[$count][0] && $lowerend < $areahash{$currchr}[$count][1]){
			#print "Low\n";
			my $low = $areahash{$currchr}[$count][1] - $lowerend;
			my $area = substr($line, 0, $low);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\t", "\t", $low, "\n";
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";	
		}
		$lowerend = $lowerend + length($line);
		$prevchar = substr($line, -1);
	}
	close GEN;

	foreach my $allchr (sort keys %areahash){
		if(!defined $areahash{$allchr}[0][2]) {next;}
		for(my $k = 0; $k < $areahash{$allchr}[0][2]; $k++){
			print OUT $allchr, "\t", $areahash{$allchr}[$k][0], "\t", $areahash{$allchr}[$k][1], "\t", $areahash{$allchr}[$k][3], "\n";
		}
	}
}
# Procedure for Masking CpG Islands
else {
	my %cgi_unsorted;
	my %cgi_sorted;
	open(CGI, "<$cpgislands") or die "cpgcounter.pl: Error: cannot open $cpgislands CpG Islands file";
	# Check for header
	my $firstline = <CGI>;		
	if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
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
	while(<CGI>){
		chomp;
		my @line = split("\t",$_);
		if($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
		# Check if duplicate
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
				$cgi_sorted{$chr}[0][2]++;
			}
		}
	}
	%cgi_unsorted = (); # Empty the unsorted hash

	# Reading input file
	$firstline = <IN>;		# Check for header
	if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
		print "No header found in $inputfile, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {} # Ignore non-standard chromosomes
		else { # Push line to hash
			$areahash{$line[0]}[0][2] = 0;
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1];
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2];
			$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0;
			#print OUT $cgi{$line[0]}[$cgi{$line[0]}[0][2]][0], "\t", $cgi{$line[0]}[$cgi{$line[0]}[0][2]][1], "\n";
			$areahash{$line[0]}[0][2]++;
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
		if(!defined $areahash{$line[0]}[0][2]) {$areahash{$line[0]}[0][2] = 0;}
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][0]=$line[1];
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][1]=$line[2];
		$areahash{$line[0]}[$areahash{$line[0]}[0][2]][3]=0;
		#print OUT $cgi{$line[0]}[$cgi{$line[0]}[0][2]][0], "\t", $cgi{$line[0]}[$cgi{$line[0]}[0][2]][1], "\n";
		$areahash{$line[0]}[0][2]++;
	}
	close IN;

	print "Scanning Genome for CpG Sites.\n";
	my $currchr;
	my $count = 0;
	my $lowerend = 0;
	my $higherend = 0;
	my $prevchar = "Z";	# Last character of previous line
	my $cpgpositioncounter = 0;	# Increments starting position for cpg island checking

	while(<GEN>){
		chomp;
		my $line = $_;
		if($line =~ /_/) {next;}
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
		if(!defined $areahash{$currchr}[0][2]) {next;}
		$line = uc($line);
		# Case 1: Line is after region, move to next region
		while($count < $areahash{$currchr}[0][2] && $lowerend > $areahash{$currchr}[$count][1]){
			$count++;
		}
		if($count == $areahash{$currchr}[0][2]) {next;}
		$higherend = $higherend + length($line);
		# Case 2: Line is before region, move to next line
		if($higherend < $areahash{$currchr}[$count][0]) {
			$lowerend = $lowerend + length($line);
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
					substr($line,$low,$high) =~ tr/CG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/CG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/CG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/CG/N/;
				}
			}
			#print "Squish\n";
			my $low = $areahash{$currchr}[$count][0] - $lowerend;
			my $high = $areahash{$currchr}[$count][1] - $areahash{$currchr}[$count][0];
			my $area = substr($line, $low, $high);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";
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
					substr($line,$low,$high) =~ tr/CG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/CG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/CG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/CG/N/;
				}
			}
			#print "All\n";
			$areahash{$currchr}[$count][3] += () = $line =~ /$cg/g;
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {
$areahash{$currchr}[$count][3]++;
print "CG\n";
}
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\n";	
			#print $line, "\t", $areahash{$currchr}[$count][3], "\n";
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
					substr($line,$low,$high) =~ tr/CG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/CG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/CG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/CG/N/;
				}
			}
			#print "High\n";
			my $high = -($higherend - $areahash{$currchr}[$count][0]);
			my $area = substr($line, $high);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			my $firstchar = substr($line, 0, 1);
			if($prevchar eq "C" && $firstchar eq "G") {
			$areahash{$currchr}[$count][3]++;
			print "CG\n";
}
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\t", "\t", $high, "\n";
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";		
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
					substr($line,$low,$high) =~ tr/CG/N/;
				}
				# Case 4: CpG Island encompasses line
				elsif($lowerend >= $cgi_sorted{$currchr}[$k][0] && $higherend <= $cgi_sorted{$currchr}[$k][1]){
					$line =~ tr/CG/N/;
				}
				# Case 5: The higher end of the line overlaps with the lower end of the CpG Island
				elsif($higherend > $cgi_sorted{$currchr}[$k][0] && $higherend < $cgi_sorted{$currchr}[$k][1]){
					my $high = -($higherend - $cgi_sorted{$currchr}[$k][0]);
					substr($line, $high) =~ tr/CG/N/;	
				}
				# Case 6: The higher end of the CpG Island overlaps with the lower end of the line
				elsif($lowerend > $cgi_sorted{$currchr}[$k][0] && $lowerend < $cgi_sorted{$currchr}[$k][1]){;
					my $low = $cgi_sorted{$currchr}[$k][1] - $lowerend;
					substr($line, 0, $low) =~ tr/CG/N/;
				}
			}
			#print "Low\n";
			my $low = $areahash{$currchr}[$count][1] - $lowerend;
			my $area = substr($line, 0, $low);
			$areahash{$currchr}[$count][3] += () = $area =~ /$cg/g;
			#print $lowerend, "\t", $higherend, "\t", $areahash{$currchr}[$count][0], "\t", $areahash{$currchr}[$count][1], "\t", "\t", $low, "\n";
			#print $area, "\t", $areahash{$currchr}[$count][3], "\n";	
		}
		$lowerend = $lowerend + length($line);
		$prevchar = substr($line, -1);
	}
	close GEN;

	foreach my $allchr (sort keys %areahash){
		if(!defined $areahash{$allchr}[0][2]) {next;}
		for(my $k = 0; $k < $areahash{$allchr}[0][2]; $k++){
			print OUT $allchr, "\t", $areahash{$allchr}[$k][0], "\t", $areahash{$allchr}[$k][1], "\t", $areahash{$allchr}[$k][3], "\n";
		}
	}
}

close OUT;

__END__
	my $cpgpositioncounter = 0;	# Starting position for each iteration

	$firstline = <IN>;		# Check for header
	if ($firstline =~ /^chr/){	# Checks to see if the first line is not a header
		print "No header found in $inputfile, processing first line.\n";
		chomp($firstline);
		my @line = split("\t",$firstline);
		if ($line[0] =~ /_/) {}
		else {
			for(my $k = $cpgpositioncounter; $k < $cgi{$line[0]}[0][2]; $k++){
				# Case 1: CpG Site is before bed array region
				# If end of in line is less than start of bed line, go to next in line
				if($cgi{$line[0]}[$k][0] > $line[2]){last;}
				# Case 2: CpG Site is after bed array region
				# If end of bed line is less than start of in line, go to next bed line
				# Increment starting position, previous lines no longer necessary for further in lines
				if($cgi{$line[0]}[$k][1] < $line[1]){
					$cpgpositioncounter++;
					next;
				}
				# Case 3: CpG Site is within bed array region
				# If in line position lies within bed line position, add methylation information if sufficient
				if($line[1] >= $cgi{$line[0]}[$k][0] && $line[2] <= $cgi{$line[0]}[$k][1]){
					$line[1] = -1;
				}
				# Case 4: CpG Site is within bed array region
				# If in line position lies within bed line position, add methylation information if sufficient
				elsif($line[1] <= $cgi{$line[0]}[$k][0] && $line[2] >= $cgi{$line[0]}[$k][1]){
					#$area_hash{$line[0]}{$line[1]}{$line[2]} = 0;
					$line[1] = $cgi{$line[0]}[$k][0]+1;
				}
				# Case 5: CpG Site is within bed array region
				# If in line position lies within bed line position, add methylation information if sufficient
				elsif($line[1] <= $cgi{$line[0]}[$k][0] && $line[2] > $cgi{$line[0]}[$k][0]){
					$line[1] = $cgi{$line[0]}[$k][0]+1;
				}
				
			}
			if($line[1] != -1) {
				#$area_hash{$line[0]}{$line[1]}{$line[2]} = 0; # Push line to hash
			}
		}
	}
	else { # If first line IS a header line
		print "Header found in $inputfile, skipping first line.\n";
	}

	while(<IN>)
	{
		chomp;
		my @line = split("\t",$_);
		if ($line[0] =~ /_/){next;}	# Ignore non-standard chromosomes
		#$area_hash{$line[0]}{$line[1]}{$line[2]} = 0;	# Push line to hash
	}
	%cgi = ();





