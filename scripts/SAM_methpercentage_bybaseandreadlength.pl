#!/usr/bin/perl
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu rgchu@ucdavis.edu
# Date: 11-19-2014
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
    1) Output File Prefix
    2,4+) Input SAM file
    3,5+) Input Sample Name (for header of output file)
" unless @ARGV > 2;

my $output = shift(@ARGV);
my $outputname = $output . ".txt";
open(OUT, ">$outputname") or die "$0 cannot open $outputname OUT outfile";
print OUT "Position";
my $cgoutput = $output . "_cgcontent.txt";
open(CG, ">$cgoutput") or die "$0 cannot open $cgoutput OUT outfile";
print CG "ReadLengths";
my %Methylation;        # Hash to hold methylation information
my %CGcontent;		# Hash to hold CG content of reads
my @inputsam;           # Input SAM Files
my @Samples;            # Input Sample Names
while(@ARGV){           # Take in the inputs
	push(@inputsam, shift(@ARGV));
	my $SampleName = shift(@ARGV);
	push(@Samples,$SampleName);
}




#############################
#  Main Loop through files  #
#############################

for(my $p = 0; $p < @Samples; $p++){    # Loop through samples
	my $SampleName = $Samples[$p];
	my $inputfile = $inputsam[$p];
	my $linecount = 0;              # Line count
	open(SAM, "<$inputfile") or die "Error: cannot open $inputfile IN infile";
	print "Processing $SampleName\n";
	while(<SAM>){ 
		chomp;
		my @line = split("\t",$_);
		if($linecount % 100000 == 0) {print "At line $linecount...\n";}
		if(defined($line[14])){} else{next;} # Skip if no data
		my $cgprocessflag = 0; # Run CG process only once
		# OUT process, methylation data
		my $methstring = substr($line[14],5);

		my $readlength = length($methstring); # Get read length

		my $string = $methstring;
		my $char = 'X'; # Methylated
		my $offset = 0;
		my $base = index($string, $char, $offset); # Find 'X' in string
		while ($base != -1) { # Until 'X' is no longer found
			# Initialize data in methylation hash
			if(!defined $Methylation{$base+1}{$SampleName}{$readlength}[0]) {
				$Methylation{$base+1}{$SampleName}{$readlength}[0] = 0;
			} 
			# Increment methylated CpG site count for this readlength at this base
			$Methylation{$base+1}{$SampleName}{$readlength}[0]++;
			$offset = $base + 1; # Increment position in string
			$base = index($string, $char, $offset); # Find next 'X'

			if($cgprocessflag == 0){
				#CG process, CG content data
				my $cgstring = substr($line[15],8,-3);
				my @num = ($cgstring =~ /[CG]/g); # Scan CGs
				my $count = @num;

				if(!defined $CGcontent{$SampleName}{$readlength}[0]) {
					$CGcontent{$SampleName}{$readlength}[0] = 0;
					$CGcontent{$SampleName}{$readlength}[1] = 0;
				}
				# Place CG count and total count into hash
				$CGcontent{$SampleName}{$readlength}[0]+=$count;
				$CGcontent{$SampleName}{$readlength}[1]+=$readlength;
				$cgprocessflag = 1; # Do not run process again for this read
			}
		}
		$string = $methstring;
		$char = 'x'; # Unmethylated
		$offset = 0;
		$base = index($string, $char, $offset); # Find 'x' in string
		while ($base != -1) { # Until 'x' is no longer found
			# Initialize data in methylation hash
			if(!defined $Methylation{$base+1}{$SampleName}{$readlength}[1]) {
				$Methylation{$base+1}{$SampleName}{$readlength}[1] = 0;
			} 
			# Increment methylated CpG site count for this readlength at this base
			$Methylation{$base+1}{$SampleName}{$readlength}[1]++;
			$offset = $base + 1; # Increment position in string
			$base = index($string, $char, $offset); # Find next 'x'

			if($cgprocessflag == 0){
				#CG process, CG content data
				my $cgstring = substr($line[15],8,-3);
				my @num = ($cgstring =~ /[CG]/g); # Scan CGs
				my $count = @num;

				if(!defined $CGcontent{$SampleName}{$readlength}[0]) {
					$CGcontent{$SampleName}{$readlength}[0] = 0;
					$CGcontent{$SampleName}{$readlength}[1] = 0;
				}
				# Place CG count and total count into hash
				$CGcontent{$SampleName}{$readlength}[0]+=$count;
				$CGcontent{$SampleName}{$readlength}[1]+=$readlength;
				$cgprocessflag = 1; # Do not run process again for this read
			}
		}
		$linecount++;
	}
	close SAM;
}



###############################
#  Loop to Print to Outfiles  #
###############################

# Create column headers with format SampleName_ReadLength_{p, c}
print "Printing to output.\n";
my %SampleOutput;
foreach my $base (sort {$a<=>$b} keys %Methylation){
	foreach my $SampleName (keys %{$Methylation{$base}}){ 
		foreach my $readlength (sort {$a<=>$b} keys %{$Methylation{$base}{$SampleName}}){
			if(!defined $SampleOutput{$SampleName}{$readlength}){
				$SampleOutput{$SampleName}{$readlength} = "\t$SampleName" . "_" . $readlength;
			}
		}
	}
}

# Column headers
# p for methylation percentage (meth/total cpg sites), c for total count of cpg sites
# Also print CG content to CG output file
foreach my $SampleName (keys %SampleOutput){ 
	foreach my $readlength (sort {$a<=>$b} keys %{$SampleOutput{$SampleName}}){
		print OUT $SampleOutput{$SampleName}{$readlength}, "_p";
		print CG $SampleOutput{$SampleName}{$readlength};
	}
}
print CG "\nCG Percentage";

foreach my $SampleName (keys %SampleOutput){ 
	foreach my $readlength (sort {$a<=>$b} keys %{$SampleOutput{$SampleName}}){
		print OUT $SampleOutput{$SampleName}{$readlength}, "_c";
		# Calculate CG percentage
		my $cgpercentage = $CGcontent{$SampleName}{$readlength}[0]/$CGcontent{$SampleName}{$readlength}[1];
		print CG "\t", $cgpercentage;
	}
}
print OUT "\n";
close CG;

# Print data, each base position at a time for all the samples
foreach my $base (sort {$a<=>$b} keys %Methylation){
	print OUT $base;
	# Print percentage methylation
	foreach my $SampleName (keys %{$Methylation{$base}}){ 
		foreach my $readlength (sort {$a<=>$b} keys %{$SampleOutput{$SampleName}}){
			# If read length less than base position, then skip
			if($base > $readlength){print OUT "\t";} 
			# Else calculate methylation information
			else{
				if(!defined $Methylation{$base}{$SampleName}{$readlength}[0]) {
					$Methylation{$base}{$SampleName}{$readlength}[0] = 0;
				} 
				if(!defined $Methylation{$base}{$SampleName}{$readlength}[1]) {
					$Methylation{$base}{$SampleName}{$readlength}[1] = 0;
				} 
				my $total = $Methylation{$base}{$SampleName}{$readlength}[0]+$Methylation{$base}{$SampleName}{$readlength}[1];
				my $PercentMeth = 0;
				# If no CpG sites, default to 0
				if($total != 0) {$PercentMeth = $Methylation{$base}{$SampleName}{$readlength}[0]/$total;}
				print OUT "\t$PercentMeth";
			}
		}
	}
	# Print total count of CpG sites
	foreach my $SampleName (keys %{$Methylation{$base}}){ 
		foreach my $readlength (sort {$a<=>$b} keys %{$SampleOutput{$SampleName}}){
			# If read length less than base position, then skip
			if($base > $readlength){print OUT "\t";}
			# Else calculate methylation information
			else{
				if(!defined $Methylation{$base}{$SampleName}{$readlength}[0]) {
					$Methylation{$base}{$SampleName}{$readlength}[0] = 0;
				} 
				if(!defined $Methylation{$base}{$SampleName}{$readlength}[1]) {
					$Methylation{$base}{$SampleName}{$readlength}[1] = 0;
				} 
				my $total = $Methylation{$base}{$SampleName}{$readlength}[0]+$Methylation{$base}{$SampleName}{$readlength}[1];
				print OUT "\t$total";
			}
		}
	}
	print OUT "\n";
}

close OUT;
