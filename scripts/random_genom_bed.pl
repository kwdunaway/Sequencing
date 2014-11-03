#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Updated: 10/31/2014
#
# This script takes a genome (ex hg 19) and outputs random positions.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "$0 needs the following arguments:
    1) output file
    2) genome
    3) number of segments
    4) length of segments
" unless @ARGV == 4;

my $outfile = shift(@ARGV);
open(OUT, ">$outfile") or die "cannot open $outfile outfile";
my $genome = shift(@ARGV);
my $segnum = shift(@ARGV);
my $length = shift(@ARGV);

my %Chroms;
my $maxpos = 0;
if($genome eq "hg19"){
	$maxpos = 249250616 - $length;
	$Chroms{"chr1"} = 249250616 - $length;
	$Chroms{"chr2"} = 243199368 - $length;
	$Chroms{"chr3"} = 198022425 - $length;
	$Chroms{"chr4"} = 191154271 - $length;
	$Chroms{"chr5"} = 180915255 - $length;
	$Chroms{"chr6"} = 171115062 - $length;
	$Chroms{"chr7"} = 159138658 - $length;
	$Chroms{"chr8"} = 146364017 - $length;
	$Chroms{"chr9"} = 141213426 - $length;
	$Chroms{"chr10"} = 135534743 - $length;
	$Chroms{"chr11"} = 135006512 - $length;
	$Chroms{"chr12"} = 133851891 - $length;
	$Chroms{"chr13"} = 115169874 - $length;
	$Chroms{"chr14"} = 107349536 - $length;
	$Chroms{"chr15"} = 102531388 - $length;
	$Chroms{"chr16"} = 90354749 - $length;
	$Chroms{"chr17"} = 81195206 - $length;
	$Chroms{"chr18"} = 78077244 - $length;
	$Chroms{"chr19"} = 59128979 - $length;
	$Chroms{"chr20"} = 63025516 - $length;
	$Chroms{"chr21"} = 48129891 - $length;
	$Chroms{"chr22"} = 51304562 - $length;
	$Chroms{"chrM"} = 16566 - $length;
	$Chroms{"chrX"} = 155270555 - $length;
	$Chroms{"chrY"} = 59373561 - $length;
}
else {die "Genome not part of built in genomes: hg19"}

my @chr;
my $chrlength = 0;
foreach my $key (sort keys %Chroms) { 
	push(@chr, $key); 
	$chrlength++;
}


##############################
# Gathering random positions #
##############################

my %BedHash;
my $count = 0;
while ($count < $segnum)
{
	my $pos = int(rand($maxpos));
	my $c = int(rand($chrlength));
	if($pos < $Chroms{$chr[$c]}){
		if(defined $BedHash{$chr[$c]}{$pos}) {next;}
		$BedHash{$chr[$c]}{$pos} = $pos + $length;
		$count++;
	}
}


############
# Printing #
############

my $number = 1;
for (my $c = 0; $c <  $chrlength; $c++){
	my $chrom = $chr[$c];
	foreach my $start (sort {$a <=> $b} values %{$BedHash{$chrom}}) {
		my $end = $start + $length;
	    print OUT $chrom , "\t", $start , "\t" , $end , "\t" , "read_" , sprintf("%07d", $number), "\n";
	    $number++;
	}
}

close OUT;
