#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: bed_compare.pl
# Version: 1.0
# Last Updated: 12-10-2012
#
# This script takes in two bed files and compares them for unique sections as well as
# overlap.  The results are printed into 3 different files.  A stats file is printed as
# well (easily loaded into Excel).
#
# Arguments:
#  (see below)
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) bed1 filename
    2) bed2 filename
    3) Bed1 Unique Name
    4) Bed2 Unique Name
    5) Bed1and2 Same Name
    6) stats filename
" unless @ARGV == 6;

my $bed1_filename = shift(@ARGV);
open(BED1IN, "<$bed1_filename") or die "cannot open $bed1_filename infile";
my $bed2_filename = shift(@ARGV);
open(BED2IN, "<$bed2_filename") or die "cannot open $bed2_filename infile";
my $bed1unique = shift(@ARGV);
my $bed1unique_filename = $bed1unique . ".bed";
open(BED1OUT, ">$bed1unique_filename") or die "cannot open $bed1unique_filename outfile";
print BED1OUT "track name=" ,$bed1unique , " description=" ,$bed1unique ;
my $bed2unique = shift(@ARGV);
my $bed2unique_filename = $bed2unique . ".bed";
open(BED2OUT, ">$bed2unique_filename") or die "cannot open $bed2unique_filename outfile";
print BED2OUT "track name=" ,$bed2unique , " description=" ,$bed2unique ;
my $bed12same = shift(@ARGV);
my $bed12same_filename = $bed12same . ".bed";
open(BED12SAMEOUT, ">$bed12same_filename") or die "cannot open $bed12same_filename outfile";
print BED12SAMEOUT "track name=" ,$bed12same , " description=" ,$bed12same ;
my $stats_filename = shift(@ARGV);
open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";
print STATS "Chromosome\t", $bed1_filename , "\t" , $bed2_filename , "\t", $bed1unique, "\t", $bed2unique, "\t", $bed12same;

my %bed1;
my %bed2;
my %Stats;



#####################################
# Load Bed files to Data Structures #
#####################################

<BED1IN>;
while (<BED1IN>){
    chomp;
    my @line = split ("\t", $_);
    push(@{$bed1{$line[0]}}, $line[1]);
    push(@{$bed1{$line[0]}}, $line[2]);
}
close BED1IN;
foreach my $chr (sort keys %bed1) {
	my $sum = 0;
	for(my $n = 0; $n < @{$bed1{$chr}}; $n++){
		$sum = $sum + $bed1{$chr}[$n + 1] - $bed1{$chr}[$n];
		$n++;
	}
   push(@{$Stats{$chr}}, $sum);	
}

<BED2IN>;
while (<BED2IN>){
    chomp;
    my @line = split ("\t", $_);
    push(@{$bed2{$line[0]}}, $line[1]);
    push(@{$bed2{$line[0]}}, $line[2]);
}
close BED2IN;

foreach my $chr (sort keys %bed2) {
	my $sum = 0;
	for(my $n = 0; $n < @{$bed2{$chr}}; $n++){
		$sum = $sum + $bed2{$chr}[$n + 1] - $bed2{$chr}[$n];
		$n++;
	}
	if(exists $Stats{$chr}){}
	else{
		push(@{$Stats{$chr}}, 0);
	}
	push(@{$Stats{$chr}}, $sum);	
	push(@{$Stats{$chr}}, 0);	
	push(@{$Stats{$chr}}, 0);	
	push(@{$Stats{$chr}}, 0);
}

foreach my $chr (sort keys %Stats) {
	if(exists $bed2{$chr}){}
	else{
		push(@{$Stats{$chr}}, 0);	
		push(@{$Stats{$chr}}, 0);	
		push(@{$Stats{$chr}}, 0);
		push(@{$Stats{$chr}}, 0);
	}
}

###########################
# Compare Data Structures #
###########################

my $bed1start = 0;
my $bed1end = 0;
my $bed2start = 99999999999;
my $bed2end = 99999999999;

foreach my $chr (sort keys %Stats) {
#	print "starting $chr \n";

	if(exists $bed1{$chr}){
		while(@{$bed1{$chr}}){
#			if($chr eq "chrY") {print $chr ,"\n" ;}
			$bed1start = shift(@{$bed1{$chr}});
			$bed1end = shift(@{$bed1{$chr}});
			if(exists $bed2{$chr}[0]){
				$bed2start = shift(@{$bed2{$chr}});
				$bed2end = shift(@{$bed2{$chr}});
			}
			else{
				my $bed2start = 99999999999;
				my $bed2end = 99999999999;
			}
		
			while($bed1start != $bed1end || $bed2start != $bed2end){
				# if both beds start same
				if($bed1start == $bed2start){
					# If both bed starts and ends are same
					if($bed1end == $bed2end){
						print BED12SAMEOUT "\n", $chr, "\t", $bed1start, "\t", $bed1end;
						$bed1start = $bed1end;
						$bed2start = $bed2end;
					}
					# If both bed starts same but bed1 ends before bed2
					elsif($bed1end < $bed2end){
						print BED12SAMEOUT "\n", $chr, "\t", $bed1start, "\t", $bed1end;
						$bed2start = $bed1end;
						if(exists $bed1{$chr}[0]){
							$bed1start = shift(@{$bed1{$chr}});
							$bed1end = shift(@{$bed1{$chr}});
						}
						else{
							$bed1start = $bed2end;
							$bed1end = $bed2end;
						}
					}
					# If both bed starts same but bed2 ends before bed1
					elsif($bed1end > $bed2end){
						print BED12SAMEOUT "\n", $chr, "\t", $bed2start, "\t", $bed2end;
						$bed1start = $bed2end;
						if(exists $bed2{$chr}[0]){
							$bed2start = shift(@{$bed2{$chr}});
							$bed2end = shift(@{$bed2{$chr}});
						}
						else{
							$bed2start = $bed1end;
							$bed2end = $bed1end;
						}
					}
					# Weirdness
					else{ die "Something went horribly wrong!!!";}
				}
				# if bed1 starts before bed2
				elsif($bed1start < $bed2start){
					# if bed1 starts and ends before bed2 starts
					if($bed1end <= $bed2start){
						print BED1OUT "\n", $chr, "\t", $bed1start, "\t", $bed1end;
						if(exists $bed1{$chr}[0]){
							$bed1start = shift(@{$bed1{$chr}});
							$bed1end = shift(@{$bed1{$chr}});
						}
						else{
							$bed1start = $bed2end;
							$bed1end = $bed2end;
						}
					}
					# if bed1 starts before bed2 starts but ends after bed2 starts
					elsif($bed1end > $bed2start){
						print BED1OUT "\n", $chr, "\t", $bed1start, "\t", $bed2start;
						$bed1start = $bed2start;
					}
					# Weirdness
					else{ die "Something went horribly wrong!!!";}
				}
				# if bed2 starts before bed1
				elsif($bed1start > $bed2start){
	
					# if bed2 starts and ends before bed1 starts
					if($bed2end <= $bed1start){
						print BED2OUT "\n", $chr, "\t", $bed2start, "\t", $bed2end;
						if(exists $bed2{$chr}[0]){
							$bed2start = shift(@{$bed2{$chr}});
							$bed2end = shift(@{$bed2{$chr}});
						}
						else{
							$bed2start = $bed1end;
							$bed2end = $bed1end;
						}
					}
					# if bed2 starts before bed1 starts but ends after bed1 starts
					elsif($bed2end > $bed1start){
						print BED2OUT "\n", $chr, "\t", $bed2start, "\t", $bed1start;
						$bed2start = $bed1start;
					}
					# Weirdness
					else{ die "Something went horribly wrong!!!";}				
				}
				# Weirdness
				else{ die "Something went horribly wrong!!!";}
			}	
		}
	}
	while(@{$bed2{$chr}}){
		$bed2start = shift(@{$bed2{$chr}});
		$bed2end = shift(@{$bed2{$chr}});
		print BED2OUT "\n", $chr, "\t", $bed2start, "\t", $bed2end;
	}
}
close BED1OUT;
close BED2OUT;
close BED12SAMEOUT;



###############################################
# Finish getting stats to print to stats file #
###############################################

open(BED1, "<$bed1unique_filename") or die "cannot open $bed1unique_filename infile";
<BED1>;
while (<BED1>){
    chomp;
    my @line = split ("\t", $_);
	$Stats{$line[0]}[2] = $Stats{$line[0]}[2] + $line[2] - $line[1];
}
close BED1;

open(BED2, "<$bed2unique_filename") or die "cannot open $bed2unique_filename infile";
<BED2>;
while (<BED2>){
    chomp;
    my @line = split ("\t", $_);
	$Stats{$line[0]}[3] = $Stats{$line[0]}[3] + $line[2] - $line[1];
}
close BED2;

open(BED12, "<$bed12same_filename") or die "cannot open $bed12same_filename infile";
<BED12>;
while (<BED12>){
    chomp;
    my @line = split ("\t", $_);
	$Stats{$line[0]}[4] = $Stats{$line[0]}[4] + $line[2] - $line[1];
}
close BED12;

my @sum = (0,0,0,0,0);
foreach my $chr (sort keys %Stats) {
	print STATS "\n" , $chr, "\t" , $Stats{$chr}[0],"\t", $Stats{$chr}[1],"\t", $Stats{$chr}[2],"\t", $Stats{$chr}[3],"\t", $Stats{$chr}[4];
	for(my $t = 0; $t < @sum; $t++){
		$sum[$t] = $sum[$t] + $Stats{$chr}[$t];
	}
} 
print STATS "\nTotal\t" , $sum[0],"\t", $sum[1],"\t", $sum[2],"\t", $sum[3],"\t", $sum[4];
close STATS;

