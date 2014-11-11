#!/usr/bin/perl
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 10-30-2014
#
# This script goes through multiple BS-Seeker2 output SAM file and returns the
# methylation statistics of all positions in the reads. This was written to determine
# if there is base position methylation bias in my data.
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
    2,3+) Input SAM file
    4,5+) Input Sample Name (for header of output file)
" unless @ARGV > 2;

my $output = shift(@ARGV);
open(OUT, ">$output") or die "$0 cannot open $output OUT outfile";
$output = $output . "60.txt";
#open(OUT60, ">$output") or die "$0 cannot open $output OUT outfile";
#my $readlength = shift(@ARGV);
print OUT "Position";
#print OUT60 "Position";
my %Methylation;
my %UnMethylation;
my @inputsam;
my @Samples;
while(@ARGV){
	push(@inputsam,shift(@ARGV));
	my $SampleName = shift(@ARGV);
	push(@Samples,$SampleName);
	print OUT "\t$SampleName" , "_Meth";
#	print OUT60 "\t$SampleName" , "_Meth";
	for(my $n = 1; $n <= 100; $n++){
		$Methylation{$SampleName}{$n}{"down"} = 0;
		$UnMethylation{$SampleName}{$n}{"down"} = 0;
		$Methylation{$SampleName}{$n}{"up"} = 0;
		$UnMethylation{$SampleName}{$n}{"up"} = 0;
	}
}




#############################
#  Main Loop through files  #
#############################

for(my $p = 0; $p < @Samples; $p++){
	my $SampleName = $Samples[$p];
	print OUT "\t$SampleName" , "_Un";
#	print OUT60 "\t$SampleName" , "_Un";
	my $inputfile = $inputsam[$p];
	my $total = 0;
	my $meth = 0;
	open(SAM, "<$inputfile") or die "Error: cannot open $inputfile IN infile";
	print "Processing $SampleName\n";
	my $source = "up";
	while(<SAM>){ 
		chomp;
		my @line = split("\t",$_);
		if(defined($line[14])){} else{next;}
		my $methstring = substr($line[14],5);

		$methstring = reverse($methstring);		
#		my $source = "up";
#		if(length($methstring) <= 60){$source = "down";}

		my $string = $methstring;
		my $char = 'X';
		my $offset = 0;
		my $result = index($string, $char, $offset);
		while ($result != -1) {
			$Methylation{$SampleName}{$result}{$source}++;
			$offset = $result + 1;
			$result = index($string, $char, $offset);
		}
		$string = $methstring;
		$char = 'x';
		$offset = 0;
		$result = index($string, $char, $offset);
		while ($result != -1) {
			$UnMethylation{$SampleName}{$result}{$source}++;
			$offset = $result + 1;
			$result = index($string, $char, $offset);
		}
	}
	close SAM;
}



###############################
#  Loop to Print to Outfiles  #
###############################

for(my $n = 1; $n <= 100; $n++){
	print OUT "\n" , $n;
#	print OUT60 "\n" , $n;
	for(my $samp = 0; $samp < @Samples; $samp++){
		print OUT "\t" , $Methylation{$Samples[$samp]}{$n}{"up"};
#		print OUT60 "\t" , $Methylation{$Samples[$samp]}{$n}{"down"};
	}
	for(my $samp = 0; $samp < @Samples; $samp++){
		print OUT "\t" , $UnMethylation{$Samples[$samp]}{$n}{"up"};
#		print OUT60 "\t" , $UnMethylation{$Samples[$samp]}{$n}{"down"};
	}
}
close OUT;
