#!/usr/bin/perl 
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Date: 5-25-2012
# Script Name: FPKM_from_GTFv2.pl
#
# Scores RPKM for the given GTF file and Bed directory/prefix
#
# Arguments:
#    1) Input BED prefix
#    2) Input GTF file
#    3) Output RPKM table file name
#    4) Total Reads (in millions)
#
################################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "useage: RPKM_from_BEDv2.pl 
    1) Input BED prefix
    2) Input GTF file
    3) Output RPKM table file name
    4) Total Reads (in millions)
" unless @ARGV == 4;

my $inputBEDprefix = shift(@ARGV);
my $GTFfilename = shift(@ARGV);
open(GTF, "<$GTFfilename") or die "cannot open $GTFfilename GTF infile";
my $OutputName = shift(@ARGV);
open(OUT, ">$OutputName") or die "cannot open $OutputName OUT outfile";
print OUT "Gene_Name\tChrom_pos\tGene_length\tStrand\tRPKM\n";
my $TotalReads = shift(@ARGV);

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
#chrom	source	feature	start	end	score	strand	frame	attribute
# ex:
#chr1    unknown exon    3204563 3207049 .       -       .       gene_id "Xkr4"; gene_name "Xkr4"; p_id "P2690"; transcript_id "NM_001011874"; tss_id "TSS1818";
# and returns the following data structure
#
# GTF{chrom}->{gene_name}->[start,stop,strand,count]
# 
# Note: count will be used later in the script but for now will be 0.

print "Loading GTFHash from GTF file.\n";
my %GTFHash;
while(<GTF>){
	chomp;
	my @line = split("\t",$_);
	my $chrom = $line[0];
	my $source = $line[1];
	my $feature = $line[2];
	my $start = $line[3];
	my $end = $line[4];
	my $score = $line[5];
	my $strand = $line[6];
	my $frame = $line[7];
	my @attribute = split("\"",$line[8]);
	my $gene_id = $attribute[1];
	my $gene_name = $attribute[3];
	my $p_id = $attribute[5];
	my $transcript_id = $attribute[7];
	my $tss_id = $attribute[9];	
	if($feature eq "exon"){
		if(exists $GTFHash{$chrom}{$gene_name}){
			if($start < $GTFHash{$chrom}{$gene_name}[0]){
				$GTFHash{$chrom}{$gene_name}[0] = $start;
			}
			if($end < $GTFHash{$chrom}{$gene_name}[0]){
				$GTFHash{$chrom}{$gene_name}[0] = $end;
			}
			if($end > $GTFHash{$chrom}{$gene_name}[1]){
				$GTFHash{$chrom}{$gene_name}[1] = $end;
			}
			if($start > $GTFHash{$chrom}{$gene_name}[1]){
				$GTFHash{$chrom}{$gene_name}[1] = $start;
			}
		}
		else{
			push(@{$GTFHash{$chrom}{$gene_name}},$start);
			push(@{$GTFHash{$chrom}{$gene_name}},$end);
			push(@{$GTFHash{$chrom}{$gene_name}},$strand);
			push(@{$GTFHash{$chrom}{$gene_name}},0);
		}
	}
}



##############################
# Calculates and prints FPKM #
##############################

while(@Chr){
	my $chr = shift(@Chr);
	my $chrom_name = "chr" . $chr;

	#opens infile
	my $inputfile = $inputBEDprefix . $chr . ".bed";
	open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";

	#Makes a StartHash for the chromosome
	print "Loading Chr$chr start positions into StartHash.\n";
	my %StartHash;  #StartHash{position}[genename1,genename2,ect]
	foreach my $gene_name (keys %{$GTFHash{$chrom_name}}) 	{ 
		push($StartHash{$GTFHash{$chrom_name}{$gene_name}[0]},$gene_name);
	}


	print "Analyzing Chromosome $chr \n";
	while(<IN>){
		chomp;
		my @InArray = split("\t",$_);
		
		my $startread = $InArray[1];  #start of read
		my $endread = $InArray[2];    #end of read

		foreach my $gene_name (keys %{$GTFHash{$chrom_name}}) 	{ 
			if($endread >= $GTFHash{$chrom_name}{$gene_name}[0] && $startread <= $GTFHash{$chrom_name}{$gene_name}[1]){
					++$GTFHash{$chrom_name}{$gene_name}[3];
			}
		}
	}
	#print GTFHash to Outfile
	foreach my $gene_name (keys %{$GTFHash{$chrom_name}}) 	{ 
		my $Chrposition = $chrom_name . ":" . $GTFHash{$chrom_name}{$gene_name}[0] . "-" . $GTFHash{$chrom_name}{$gene_name}[1];
		my $genelength = $GTFHash{$chrom_name}{$gene_name}[1] - $GTFHash{$chrom_name}{$gene_name}[0];
		my $RPKMscore = $GTFHash{$chrom_name}{$gene_name}[3] / $genelength;
		$RPKMscore = $RPKMscore / $TotalReads;
		$RPKMscore = sprintf("%.4f", $RPKMscore);		
		print OUT "$gene_name\t$Chrposition\t$genelength\t$GTFHash{$chrom_name}{$gene_name}[2]\t$RPKMscore\n";
	}
}
