#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-25-2012
# Script Name: FPKM_from_GTF.pl
#
# Scores FPKM for the given GTF file
#
# Arguments:
#    1) Input FPKMWIG prefix
#    2) Input GTF file
#    3) Output FPKM file 
#    4) Column name for table
#
###############################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: FPKM_from_GTF.pl 
    1) Input FPKMWIG prefix
    2) Input GTF file
    3) Output FPKM summary table file name 
    4) Column name for table
" unless @ARGV == 4;

my $inputFPKMprefix = shift(@ARGV);
my $GTFfilename = shift(@ARGV);
open(GTF, "<$GTFfilename") or die "cannot open $GTFfilename GTF infile";
my $OutputName = shift(@ARGV);
open(OUT, ">$OutputName") or die "cannot open $OutputName OUT outfile";
my $ColumnName = shift(@ARGV);

print OUT "Feature_Name\tFull_position\tSub_Positions\tLength\tStrand\t$ColumnName\n";


my @Chr;             # array that contains all the the names of the mouse chromosomes
for (my $n = 1; $n< 20; $n++){
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");



########################
# GTF structure create #
########################

# Takes a GTF file with lines like this:
#chrom	source	feature	start	end	score	strand
# ex:
#chr1    CpG_Island CpG_100    3204563 3207049 .       +
# and returns the following data structure
#
# GTF{chrom}->{feature}->{number}->[start,stop,length,strand,summed FPKM over length]
# 
# Note: FPKM will be used later in the script but for now will be 0.

my $source;

print "Loading GTFHash from GTF file.\n";
my %GTFHash;
while(<GTF>){
	chomp;
	my @line = split("\t",$_);
	my $chrom = $line[0];
	$source = $line[1];
	my $feature = $line[2];
	my $start = $line[3];
	my $end = $line[4];
	my $score = $line[5];
	my $strand = $line[6];
	
	# Corrects if start is greater than end
	if($start > $end){
		my $temp = $end;
		$end = $start;
		$start = $temp;
	}
	my $length = $end - $start;
	my $flag = 1;
	while($flag > 0){
		if(exists $GTFHash{$chrom}{$feature}{$flag}){
			$flag++;
		}
		else{
			push(@{$GTFHash{$chrom}{$feature}{$flag}},$start);
			push(@{$GTFHash{$chrom}{$feature}{$flag}},$end);
			push(@{$GTFHash{$chrom}{$feature}{$flag}},$length);
			push(@{$GTFHash{$chrom}{$feature}{$flag}},$strand);
			push(@{$GTFHash{$chrom}{$feature}{$flag}},0);
			$flag = 0;
		}
	}
}
close GTF;



##############################
# Calculates and prints FPKM #
##############################

while(@Chr){
#	my %FPKMHash;  #holds FPKM from IN file for small selected area
	my $chr = shift(@Chr);

	#opens infile and removes first two lines from it (header lines)
	my $inputfile = $inputFPKMprefix . $chr . ".wig";
	open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";
	<IN>; <IN>;

	my %FeatureHash;
	my $chrom_name = "chr" . $chr;

	print "Loading Chr$chr start positions into FeatureHash.\n";
	my %StartHash;
	#Makes a StartHash and a FeatureHash for the chromosome
	foreach my $feature (keys %{$GTFHash{$chrom_name}}) 	{ 
		my $feature_start = 999999999999;
		my $feature_end = 0;
		foreach my $segnumber (keys %{$GTFHash{$chrom_name}{$feature}}) 	{ 
			if($GTFHash{$chrom_name}{$feature}{$segnumber}[0] < $feature_start){
				$feature_start = $GTFHash{$chrom_name}{$feature}{$segnumber}[0];
			}
			if($GTFHash{$chrom_name}{$feature}{$segnumber}[1] > $feature_end){
				$feature_end = $GTFHash{$chrom_name}{$feature}{$segnumber}[1];
			}		
		}
		#Adds start and end to Feature Hash
		push(@{$FeatureHash{$feature}},$feature_start);
		push(@{$FeatureHash{$feature}},$feature_end);
		
		#Adds start position and corresponding 
		#Note, a single starting position can have multiple features associated with it
		push(@{$StartHash{$feature_start}},$feature);
	}

	my @SortedFeatureArray;
	# Fills @SortedFeatureArray with list (in order of starting position)
	foreach my $startpos (sort { $a <=> $b } keys %StartHash) 	{
		while(exists $StartHash{$startpos}[0]){
			push(@SortedFeatureArray,$StartHash{$startpos}[0]);
			shift(@{$StartHash{$startpos}});
		}
	}
	
	undef %StartHash;

	my $lastendpos = 0;

	print "Analyzing Chromosome $chr \n";
	while(<IN>){
		chomp;
		my @InArray = split("\t",$_);
		my $currentpos = $InArray[0];
		my $currentheight = $InArray[1];

		# Ends the loop if there are no more features in this array
		if(exists $SortedFeatureArray[0]){}
		else {last;}
		
		my $sortedfeaturepointer = 0;
		# while current position is not before the start of the current hash
		while($FeatureHash{$SortedFeatureArray[$sortedfeaturepointer]}[0] < $currentpos){
			my $feature = $SortedFeatureArray[$sortedfeaturepointer];
			
			#If current position is past the end of the feature, removes feature from SortedFeatureArray
			if($currentpos > $FeatureHash{$feature}[1]){
				shift(@SortedFeatureArray);
#				my $feature_finished = shift(@SortedFeatureArray);
#				print "Finished Feature: " , $feature_finished , "\t", $FeatureHash{$feature}[1] , "\n";
				if(exists $SortedFeatureArray[0]){next;}
				else {last;}
			}		
			
			#If current position is before the feature start, exit loop
			if($currentpos < $FeatureHash{$feature}[0]){
				last;
			}
						
			#Go through feature area, adding FPKM to areas that overlap with it
			my $segnumber = 1;
			while(exists $GTFHash{$chrom_name}{$feature}{$segnumber}){
				if($currentpos > $GTFHash{$chrom_name}{$feature}{$segnumber}[0] && $currentpos < $GTFHash{$chrom_name}{$feature}{$segnumber}[1]){
					$GTFHash{$chrom_name}{$feature}{$segnumber}[4] = $GTFHash{$chrom_name}{$feature}{$segnumber}[4] + $currentheight;
				}
				++$segnumber;
			}			
			++$sortedfeaturepointer;
		}
	}
	close IN;
	
	print "Printing $chrom_name to Outfile\n\n";
	#Print FPKM's to Outfile
	foreach my $feature (sort keys %{$GTFHash{$chrom_name}}) 	{ 
		my $featurelength = 0;
		my $FPKM = 0;
		my $Chr_total_position = $chrom_name . ":" . $FeatureHash{$feature}[0] . "-" . $FeatureHash{$feature}[1];
		my $Chrpositions = "";
		my $loopcount = 0;
		foreach my $segnumber (keys %{$GTFHash{$chrom_name}{$feature}}) {
			if ($loopcount > 0){
				$Chrpositions = $Chrpositions . " / ";
			}
			++$loopcount;
			$featurelength = $featurelength + $GTFHash{$chrom_name}{$feature}{$segnumber}[2];
			$FPKM = $FPKM + $GTFHash{$chrom_name}{$feature}{$segnumber}[4];
			$Chrpositions = $chrom_name . ":" . $GTFHash{$chrom_name}{$feature}{$segnumber}[0] . "-" . $GTFHash{$chrom_name}{$feature}{$segnumber}[1];
		}
		if($featurelength > 0){
			$FPKM = (1000 * $FPKM) / $featurelength;
			$FPKM = sprintf("%.4f", $FPKM);	
		} 
		else{
			$FPKM = 0;
		}
		# Print Feature	totalposition	positions	Feature Length(s)	strand	FPKM
		print OUT "$feature\t$Chr_total_position\t$Chrpositions\t$featurelength\t$GTFHash{$chrom_name}{$feature}{1}[3]\t$FPKM\n";
	}	
}
close OUT;
