#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Update Date: 6-3-2014
# Version: 1.0
#
# Takes SAM output from BS_Seeker2 and creates percentage methylation BED files that
# can be uploaded to the UCSC genome browser or further analyzed through StochHMM.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following parameters:
    1) Input sorted SAM file
    2) Output files prefix (folder and prefix)
    3) Bed track prefix
    4) UCSC genome version (ex: hg19)
    5) Methylation type (CG, CHG, or CHH)
    6) Strand (combined, positive, or negative)
" unless @ARGV == 6;

my $infile = shift(@ARGV);
open(IN, "<$infile") or die "cannot open $infile infile";
my $outprefix = shift(@ARGV);
my $bedprefix = shift(@ARGV);
my $genome_version = shift(@ARGV);
my $meth_type = shift(@ARGV);
unless (($meth_type eq "CG") || ($meth_type eq "CHG") || ($meth_type eq "CHH")) {
    die "Methylation type $meth_type is not one of:  CG, CHG, or CHH\n\n";
}
my $strand_type = shift(@ARGV);
unless (($strand_type eq "combined") || ($strand_type eq "positive") || ($strand_type eq "negative")) {
    die "Strand type $strand_type is not one of: combined, positive, or negative\n\n";
}

# Global Variables
my $currentchrom = "Not Set Yet";
my %Methylation;
my $prevstart = 0;
my $prevstrand = "+";

#Columns for formatting SAM files
my $chrc = 2;
my $startc = 3;
my $methc = 14;
my $strandc = 11;



#############
# Main Loop #
#############

while(<IN>){
	my @line = split("\t", $_);
	my $chrom = $line[$chrc];
	if ($chrom =~ /_/) {next;}
	my $start = $line[$startc];
	my $methstring = substr $line[$methc], 5;
	my $strand = substr $line[$strandc], 5,1;

	# If duplicate line, skip
	if($prevstart == $start && $prevstrand eq $strand) {next;}
	
	# If only looking at positive strand
	if($strand eq "-" && $strand_type eq "positive") {next;}
	# If only looking at negative strand
	if($strand eq "+" && $strand_type eq "negative") {next;}

	# On next chromosome
	if($chrom ne $currentchrom){
		print "Starting " , $chrom , "\n";
		if($currentchrom ne "Not Set Yet"){
			my $outfile = $outprefix . $currentchrom . ".bed";
			open(OUT, ">$outfile") or die "cannot open $outfile outfile";
			print OUT "track name=" , $bedprefix, $currentchrom, " description=" , $bedprefix, $currentchrom, " useScore=0 itemRgb=On db=" , $genome_version , "\n";
			foreach my $posstart (sort { $a <=> $b } keys %Methylation) {
				#Example print format
				#chr10   51332   51333   0.50-2  0       +       0       0       27,74,210
				my $posend = $posstart + 1;
				my $methperc = 0;
				my @methraw = split("",$Methylation{$posstart});
				my $methnum = @methraw;
				while(@methraw){
					$methperc += $methraw[0];
					shift(@methraw);
				}
				$methperc = $methperc / $methnum;
				$methperc = sprintf("%.2f", $methperc);
				my $color = "0,0,0"; #black
				if ($methperc > 0 && $methperc <= .6) {$color = "27,74,210";} #blue
				elsif ($methperc > .6 && $methperc <= .8) {$color = "27,210,57";} #green
				elsif ($methperc > .8) {$color = "210,27,27";} #red
				print OUT $currentchrom , "\t" , $posstart , "\t" , $posend , "\t" , $methperc , "-", $methnum , "\t" , "0\t+\t0\t0\t" , $color , "\n";
			}
		}
		close OUT;
		$currentchrom = $chrom;
		%Methylation = ();
	}
	

	# Find X
	my $offset = 0;
	my $position = 0;
	while ($position >= 0)
	{
  		$position = index($methstring, 'X', $offset);
  		if($position == -1) {last;}
 		if($strand eq "+"){
			my $startpos = $start + $position;
			if(defined $Methylation{$startpos}) {$Methylation{$startpos} = $Methylation{$startpos} . "1";}
			else {$Methylation{$startpos} = "1";}
 		}
 		elsif($strand eq "-"){
			my $startpos = $start + length($methstring) - $position -2;
			if(defined $Methylation{$startpos}) {$Methylation{$startpos} = $Methylation{$startpos} . "1";}
			else {$Methylation{$startpos} = "1";}
 		}
 		else { die "Strand not + or - but $strand \n";}
		$offset=$position+1;
	}

	# Find x
	$offset = 0;
	$position = 0;
	while ($position >= 0)
	{
  		$position = index($methstring, 'x', $offset);
  		if($position == -1) {last;}
 		if($strand eq "+"){
			my $startpos = $start + $position;
			if(defined $Methylation{$startpos}) {$Methylation{$startpos} = $Methylation{$startpos} . "0";}
			else {$Methylation{$startpos} = "0";}
 		}
 		elsif($strand eq "-"){
			my $startpos = $start + length($methstring) - $position -2;
			if(defined $Methylation{$startpos}) {$Methylation{$startpos} = $Methylation{$startpos} . "0";}
			else {$Methylation{$startpos} = "0";}
 		}
 		else { die "Strand not + or - but $strand \n";}
		$offset=$position+1;
	}
	$prevstart = $start;
	$prevstrand = $strand;
}

my $outfile = $outprefix . $currentchrom . ".bed";
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
print OUT "track name=" , $bedprefix, $currentchrom, " description=" , $bedprefix, $currentchrom, " useScore=0 itemRgb=On db=" , $genome_version , "\n";
foreach my $posstart (sort { $a <=> $b }  keys %Methylation) {
	#Example print format
	#chr10   51332   51333   0.50-2  0       +       0       0       27,74,210
	my $posend = $posstart + 1;
	my $methperc = 0;
	my @methraw = split("",$Methylation{$posstart});
	my $methnum = @methraw;
	while(@methraw){
		$methperc += $methraw[0];
		shift(@methraw);
	}
	$methperc = $methperc / $methnum;
	$methperc = sprintf("%.2f", $methperc);
	my $color = "0,0,0"; #black
	if ($methperc > 0 && $methperc <= .6) {$color = "27,74,210";} #blue
	elsif ($methperc > .6 && $methperc <= .8) {$color = "27,210,57";} #green
	elsif ($methperc > .8) {$color = "210,27,27";} #red
	print OUT $currentchrom , "\t" , $posstart , "\t" , $posend , "\t" , $methperc , "-", $methnum , "\t" , "0\t+\t0\t0\t" , $color , "\n";
}
close OUT;

#create subroutines
#print chromosome
#check line for methylation and add to hash (make strand and methy type dependent)
