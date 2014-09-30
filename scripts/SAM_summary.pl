#!/usr/bin/perl
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 9-23-2014
#
# This script goes through multiple BS-Seeker2 output SAM file and returns the
# methylation statistics of all reads.
#
# Output:
# TotalCpGsInRead MethylatedCpGs Sample1counts Sample2counts...
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following parameters:
    1) Output File Name
    2,4+) Input SAM file
    3,5+) Input Sample Name (for header of output file)
" unless @ARGV > 2;

my $output = shift(@ARGV);
open(OUT, ">$output") or die "$0 cannot open $output OUT outfile";
print OUT "CpGsInRead", "\t" , "MethCpGs";
my @inputsam;
my @Samples;
while(@ARGV){
	push(@inputsam,shift(@ARGV));
	push(@Samples,shift(@ARGV));
}

#Format: $StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName} = count;
my %StrandMethylation;

#############################
#  Main Loop through files  #
#############################

for(my $p = 0; $p < @Samples; $p++){
	my $SampleName = $Samples[$p];
	print OUT "\t$SampleName";
	my $inputfile = $inputsam[$p];
	my $total = 0;
	my $meth = 0;
	open(SAM, "<$inputfile") or die "Error: cannot open $inputfile IN infile";
	while(<SAM>){ 
		chomp;
		my @line = split("\t",$_);
		my $methstring = substr($line[14],5);
		my $Methylated_CpGs = ($methstring =~ tr/X//);
		my $Total_CpGs = ($methstring =~ tr/x//) + $Methylated_CpGs;
		if ($Total_CpGs == 0){next;}
		$total += $Total_CpGs;
		$meth += $Methylated_CpGs;
		if(defined $StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName}){
			$StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName}++;
		}
		else{$StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName} = 1;}
	}
	close SAM;
	my $avgmeth = $meth / $total;
	print "Average methylation of all reads for sample $SampleName\t$avgmeth\n";
}



###############################
#  Loop to Print to Out File  #
###############################

foreach my $Total_CpGs (sort { $a <=> $b } keys (%StrandMethylation) ){
	foreach my $Methylated_CpGs (sort { $a <=> $b } keys (%{$StrandMethylation{$Total_CpGs}}) ){
		print OUT "\n" , $Total_CpGs , "\t" , $Methylated_CpGs;
		for(my $p = 0; $p < @Samples; $p++){
			my $SampleName = $Samples[$p];
			if(defined $StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName}){
	    		print OUT "\t" , $StrandMethylation{$Total_CpGs}{$Methylated_CpGs}{$SampleName};
	    	}
	    	else {print OUT "\t0";}
    	}
   	}
}
close OUT;
