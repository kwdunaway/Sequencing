#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 10-3-2014
#
# This script will yield a summary of every CpG in the PercentMethfile.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) Output file
    2,4+) Input Percent Methylation Folder Prefix (no chr)
    3,5+) Input Sample Name (for header of output file)
" unless @ARGV > 2;

my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
my @filenames;			# Input Percent Methylation Folder Prefix (Ex: /home/user/Permeth/Permeth_)
my @headernames;		# Input Percent Methylation Header for Column (Ex: Sample1)
while(@ARGV){
	push(@filenames,shift(@ARGV));
	push(@headernames,shift(@ARGV));
}

my @Master_chrom;
for (my $n = 1; $n < 23; $n++){
	push(@Master_chrom, $n);
}
push(@Master_chrom, "X");
push(@Master_chrom, "Y");
my $alltotalCpGs = 0;
my $alltotalassays = 0;
my $alltotalCs = 0;
my $allpercentSum = 0;

# Format of datastructure: $MethHash{SampleName}{chr}{key} = value;
my %MethHash;


#############
# Main Loop #
#############

for(my $i = 0; $i < @filenames; $i++){
	my $SampleName = $headernames[$i];
	$MethHash{$SampleName}{"all"}{"totalCpGs"} = 0;
	$MethHash{$SampleName}{"all"}{"totalassays"} = 0;
	$MethHash{$SampleName}{"all"}{"totalCs"} = 0;
	$MethHash{$SampleName}{"all"}{"percentSum"} = 0;

	my @chrom = @Master_chrom;
	while(@chrom){
		my $chr = shift(@chrom);
		my $chromosome = "chr" . $chr;
		my $filename = $filenames[$i] . $chr . ".bed";	# Create File Name
		open (IN, "<$filename") or die "$0: Error: Couldn't open chromosome file $filename\n";
		print "Analyzing $filename\n";

		my $totalCpGs = 0;
		my $totalassays = 0;
		my $totalCs = 0;
		my $percentSum = 0;
		#remove header
		<IN>;

		while (<IN>){
		    chomp;
		    my @line = split ("\t", $_);
		    my @tmp = split ("-", $line[3]);
		    my $methperc = $tmp[0];
		    my $coverage = $tmp[1];
		    my $Cs = int(($methperc * $coverage) +0.5);
			$totalCpGs++;
			$totalassays+= $coverage;
			$totalCs+= $Cs;
			$percentSum += $methperc;
		}
		close IN;
		$MethHash{$SampleName}{$chromosome}{"totalCpGs"} = $totalCpGs;
		$MethHash{$SampleName}{$chromosome}{"totalassays"} = $totalassays;
		$MethHash{$SampleName}{$chromosome}{"totalCs"} = $totalCs;
		$MethHash{$SampleName}{$chromosome}{"percentSum"} = $percentSum;
		$MethHash{$SampleName}{$chromosome}{"avgmeth"} = $totalCs/$totalassays;
		$MethHash{$SampleName}{$chromosome}{"avgcoverage"} = $totalassays/$totalCpGs;
		$MethHash{$SampleName}{$chromosome}{"percentAvg"} = $percentSum/$totalCpGs;
		
		$MethHash{$SampleName}{"all"}{"totalCpGs"} += $totalCpGs;
		$MethHash{$SampleName}{"all"}{"totalassays"} += $totalassays;
		$MethHash{$SampleName}{"all"}{"totalCs"} += $totalCs;
		$MethHash{$SampleName}{"all"}{"percentSum"} += $percentSum;

	}
	$MethHash{$SampleName}{"all"}{"avgmeth"} = $MethHash{$SampleName}{"all"}{"totalCs"} / $MethHash{$SampleName}{"all"}{"totalassays"};
	$MethHash{$SampleName}{"all"}{"avgcoverage"} = $MethHash{$SampleName}{"all"}{"totalassays"} / $MethHash{$SampleName}{"all"}{"totalCpGs"};
	$MethHash{$SampleName}{"all"}{"percentAvg"} = $MethHash{$SampleName}{"all"}{"percentSum"} / $MethHash{$SampleName}{"all"}{"totalCpGs"};
}


#my @Variables = ("totalCpGs","totalassays","totalCs","percentSum","avgmeth","avgcoverage","percentAvg");

foreach my $variable (keys (%{$MethHash{$headernames[0]}{"all"}})) {
	print OUT $variable , "\n";
	print OUT "Chromosome";
	for(my $i = 0; $i < @filenames; $i++){
		my $SampleName = $headernames[$i];
		print OUT "\t" , $SampleName;
	}
	print OUT "\n";
	for(my $c = 0; $c < @Master_chrom; $c++){
		my $chromosome = "chr" . $Master_chrom[$c];
		print OUT $chromosome;
		for(my $i = 0; $i < @filenames; $i++){
			my $SampleName = $headernames[$i];
			print OUT "\t" , $MethHash{$SampleName}{$chromosome}{$variable};
		}
		print OUT "\n";
	}
	for(my $i = 0; $i < @filenames; $i++){
		my $SampleName = $headernames[$i];
		print OUT "\t" , $MethHash{$SampleName}{"all"}{$variable};
	}
	print OUT "\n\n";
}
close OUT;


__END__
print OUT "\nTotal CpGs";
for(my $i = 0; $i < @filenames; $i++){
	my $SampleName = $headernames[$i];
	print OUT "\t" , $MethHash{$SampleName}{"all"}{"totalCs"};
}
print OUT "\nTotal AvgCoverage";
for(my $i = 0; $i < @filenames; $i++){
	my $SampleName = $headernames[$i];
	print OUT "\t" , $MethHash{$SampleName}{"all"}{"avgcoverage"};
}
print OUT "\nTotal AvgMeth";
for(my $i = 0; $i < @filenames; $i++){
	my $SampleName = $headernames[$i];
	print OUT "\t" , $MethHash{$SampleName}{"all"}{"avgmeth"};
}
print OUT "\nTotal PercentAvg";
for(my $i = 0; $i < @filenames; $i++){
	my $SampleName = $headernames[$i];
	print OUT "\t" , $MethHash{$SampleName}{"all"}{"percentAvg"};
}

my @chrom = @Master_chrom;

while(@chrom){
	my $chr = shift(@chrom);
	my $chromosome = "chr" . $chr;
	print OUT $chromosome;
	for(my $i = 0; $i < @filenames; $i++){
		my $SampleName = $headernames[$i];
		print OUT "\t" , $MethHash{$SampleName}{$chromosome}{"percentAvg"};
	}
	print OUT "\n";
}
close OUT;

__END__
my $allavgmeth = $alltotalCs/$alltotalassays;
my $allavgcoverage = $alltotalassays/$alltotalCpGs;
print "Total for Experiment:
	CpGs:        \t$alltotalCpGs
	Total Assays:\t$alltotalassays
	Total C's:   \t$alltotalCs
	Avg Meth:    \t$allavgmeth
	Avg Coverage \t$allavgcoverage\n\n";


