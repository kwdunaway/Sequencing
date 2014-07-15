#!/usr/bin/perl
BEGIN {push @INC, "/data/scratch/programs/perl_script";}
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 7-11-2014
# Script Name: strandbias.pl
#
# This script scans CpG's of reads in windows of the input file
#
# Arguments: See Below
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "strandbias.pl needs the following parameters:
    1) Minimum Threshold for Number of CpG's
    2) Input BED File or type 'all' for entire genome (10kbp increments)
    3) CpG Islands BED File for masking or type 'nomask'
    4) Output File Name
    5+) Input SAM File(s)
    IMPORTANT: Include the full path name to the input files
" unless @ARGV >= 5;

my $threshold = shift(@ARGV);
my $inputbed = shift(@ARGV);
my $cpgislands = shift(@ARGV);
my $output = shift(@ARGV);
my @inputsam = @ARGV;

open(OUT, ">$output") or die "Error: strandbias.pl: cannot open $output OUT outfile";

print "\nStrand Bias Script\n\n";

#########################
#     Creating Hash     #
#########################

my @methinfo; # 2-d Array storing the output
my $linecount = 0; # Number of lines for the output
my $samples = $#inputsam + 1;

if($inputbed eq "all")
{
	print "all\n";
}
else
{
	open(IN, "<$inputbed") or die "Error: strandbias.pl: cannot open $inputbed IN infile";
	while(<IN>){ 
		chomp;
		my @line = split("\t",$_);
		$methinfo[$linecount][0][0][0] = $line[0];	#chromosome
		$methinfo[$linecount][0][0][1] = $line[1];	#start
		$methinfo[$linecount][0][0][2] = $line[2];	#end
		for(my $i = 0; $i < $samples; $i++){		# for every sample
			$methinfo[$linecount][$i][0][0] = 0;	# initialize #meth count
			# [line count] [SAM] [read] [0:#CpG 1-#:#methylated]
		}
		$linecount++;
	}

	close IN;
}

##############################################
#     Processing SAM Files and Filtering     #
##############################################

if($cpgislands eq "nomask")
{
	for(my $i = 0; $i < $samples; $i++){		# for every sample
		my $sam = shift(@inputsam);
		print $sam, "\n";
		open(SAM, "<$sam") or die "Error: strandbias.pl: cannot open $sam SAM infile";
		while(<SAM>)
		{
			
		}

		close SAM;
	}
}
else	#CpG Island BED File was input, mask CpG Island-containing reads
{
	open(CPG, "<$cpgislands") or die "Error: strandbias.pl: cannot open $cpgislands CPG infile";

	for(my $i = 0; $i < $samples; $i++){		# for every sample
		my $sam = shift(@inputsam);
		print $sam, "\n";
		open(SAM, "<$sam") or die "Error: strandbias.pl: cannot open $sam SAM infile";
		while(<SAM>)
		{
		
		}

		close SAM;
	}

	close CPG;
}

close OUT;

__END__



