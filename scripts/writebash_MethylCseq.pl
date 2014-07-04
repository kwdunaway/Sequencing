#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: 
# Version: 0.1
# Last Updated: 
#
# <brief description of script's function>
#
# Arguments:
#   <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Name of outbash file you want to create
    2,4+) Name of sample (to be used as multiple file naming)
    3,5+) Path to folder with all gzip fastq files (without last /) (uses all *.gz files in folder)
" unless @ARGV > 1 && @ARGV %2 == 1;

my $output_filename = shift(@ARGV);
open(OUT, ">$output_filename") or die "cannot open $output_filename outfile";

#my $SampName = "JLKD009";
#my $FASTQ_folder = "Sample_JLKD009";

#Must be changed if aligning to genome other than hg18 -e=80 -m=2
#Path to BS_Seeker compiled Reference Genome
my $BSRefGenome = "/data/scratch/genomes/hg18/reference_genome/";
my $genome = "hg18";
my $e = "80";
my $m = "2";

#Must change depending on computer
#Path to BS_Seeker
my $BSPath = "/data/scratch/programs/tuxedo/BS_Seeker/";
#Path to Bowtie
my $BowtiePath = "/data/scratch/programs/tuxedo/bowtie-0.12.7/";
#Path to Samtools
my $SamtoolsPath = "/data/scratch/programs/tuxedo/samtools/";
#Path (complete) to out2browserview (makes Percent Methylation files
my $OuttobrowserviewPath = "/data/scratch/programs/perl_script/BS_Seeker_out2browserview_v4.pl";

###################
# Writing to BASH #
###################

#set PATH for program calls
print OUT "PATH=\$PATH:" , $BowtiePath , ":" , $SamtoolsPath , "\n\n";

while(@ARGV){
	my $SampName = shift(@ARGV);
	my $FASTQ_folder = shift(@ARGV);
	
	#prints Sample information
	print OUT "#SampleName=$SampName FASTQ_folder=$FASTQ_folder\n";

	#unzip and combine fastq files
	print OUT "gunzip -c " , $FASTQ_folder  , "/*.gz | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' |   grep -v \"^--\$\" > " , $FASTQ_folder , "/" , $SampName , ".fq\n";

	#split fastq files into 20 mil chunks
	print OUT "split -l 80000000 " , $FASTQ_folder , "/" , $SampName , ".fq " , $FASTQ_folder , "/" , $SampName , "_\n";

	#align reads using -e 80 -m 2
	for (my $n = 'aa'; $n ne "ak"; ++$n){
		print OUT "python " , $BSPath , "BS_Seeker.py -t N -e " , $e , " -m " , $m , " -a " , $BSPath , "adapter2.txt -d " , $BSRefGenome , " -p " , $BowtiePath , " -i " , $FASTQ_folder , "/" , $SampName , "_" , $n , " -o " , $SampName , "_" , $n , ".txt\n";
	}

	#combine aligned reads to one file
	my $BSoutfilename = $SampName . "_BSOUT_e" . $e . "m" .$m . ".txt";
	print OUT "cat " , $SampName , "_a* > " , $BSoutfilename , "\n";
	print OUT "cat log_" , $SampName , "_a* > log_" , $BSoutfilename , "\n";

	#create permeth files
	print OUT "perl " , $OuttobrowserviewPath , " " , $BSoutfilename , " Sorted_" , $SampName , " " , $SampName , " " , $genome , " n\n\n";
}