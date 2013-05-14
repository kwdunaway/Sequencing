#!/usr/bin/perl
BEGIN {push @INC, "/home/kwdunaway/perl_script";}
use strict; use warnings;
use SeqProcess;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 5-12-2013
# Script Name: ChIP_SeqProcess.pl
#
# This script processes the subroutines in SeqProcess.pm.
#
# Arguments: See Below
#
##########################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "ChIP_SeqProcess.pl needs the following parameters:
    1) Folder path to contain data (e.g. /home/user/Folder/)(No spaces in path)
    2) Raw data file folder (the extension is .fq.gz (Fastq gzipped) and directory only contains .fq.gz) 
    3) File prefix for created files (general name for new files)
    4) Species (mm9 or hg18)
    5) Final Read length (for reads to be extended to)
    6) WIG Track Color (in RRR,GGG,BBB format)
    7) Maximum Duplicate Reads (1 for no duplicate reads)
" unless @ARGV == 7;

my $ExperimentTopDir = shift(@ARGV);
my $rawfqfolder = shift(@ARGV);
my $FilePrefix = shift(@ARGV);
my $Species = shift(@ARGV);
my $FinalReadLength = shift(@ARGV);
my $WIGTrackColor = shift(@ARGV);
my $MaxDupReads = shift(@ARGV);

# Error Checking
die "Error: 6) Maximum duplicate reads should be at least 1.\n" unless $MaxDupReads > 0;

# Adjusting Folder Names
$ExperimentTopDir = $ExperimentTopDir . "\/" unless $ExperimentTopDir =~ m/\/$/; # Ensures ending "/"
$rawfqfolder = $rawfqfolder . "\/" unless $rawfqfolder =~ m/\/$/; # Ensures ending "/"

my $commandline = ""; #inputs for command line


###########################
# Computer specific paths #
###########################

# Bowtie path
my $addtoPATH = ":/home/kwdunaway/tuxedo/bowtie-0.12.7";
# Tophat path
$addtoPATH = $addtoPATH . ":/home/kwdunaway/tuxedo/tophat-1.4.1.Linux_x86_64";
# Samtools path
$addtoPATH = $addtoPATH . ":/home/kwdunaway/tuxedo/samtools";
# Cufflinks path
$addtoPATH = $addtoPATH . ":/home/kwdunaway/tuxedo/cufflinks-1.3.0.Linux_x86_64";

my $refpath = "";

if($Species eq "mm9"){
	$refpath = "/work/kwdunaway/mm9/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome";
} elsif ($Species eq "hg18"){
	$refpath = "/work/kwdunaway/hg18/hg18";
} else {die "Species needs to be either mm9 or hg18, not: " , $Species;}

####################
# Global Variables #
####################

my $filtered_fastq = $rawfqfolder . "filtered.fq";
my $fastqgzfile = $rawfqfolder . "*.gz";
my $ReadLength = SeqProcess::fastqgz_readlength($fastqgzfile);  # Determines Read Length from zipped raw files
my $nonalignedreadsfile = $ExperimentTopDir . $FilePrefix . "_NonAligned.fq";
my $alignedpreseparationfile = $ExperimentTopDir . $FilePrefix . "_alignedpreseparation.txt";
my $uniqalignedreadsfile = $ExperimentTopDir . $FilePrefix . "_Uniq.txt";
my $repalignedreadsfile = $ExperimentTopDir . $FilePrefix . "_Repeat.txt";
my $bedtowigfiles = $ExperimentTopDir . $FilePrefix . "_bed/" . $FilePrefix;
my $fpkmwigfiles =  $ExperimentTopDir . $FilePrefix . "_FPKMWIG/" . $FilePrefix . "_FPKM";
my $visfpkmwig = $ExperimentTopDir . $FilePrefix . "_VisFPKMWIG/" . $FilePrefix . "_VisFPKMWIG";


###########################################################################
#                          Bowtie and Pre-bowtie                          #
#                                                                         #
# (1) Pathing                                                             #
# (2) Prepare temp files for bowtie                                       #
# (3) Bowtie (Separates reads between aligned and non-aligned)            #
# (4) Remove temp files                                                   #
###########################################################################

# (1) Header and path modifier
SeqProcess::add_path($addtoPATH);

# (2) Unzip zipped files, filter them, and combine into one .fq file
#SeqProcess::filter_zip($rawfqfolder);

# (3) Make folder for experiment and run Bowtie (separates reads between aligned and non-aligned)
SeqProcess::run_bowtie($ExperimentTopDir, $FilePrefix, $refpath, $filtered_fastq);

# (4) Remove made files
#print "Removing $filtered_fastq\n";
#`rm $filtered_fastq`;


########################################################################################
#                                Post-bowtie                                           #
#                                                                                      #
# (5) Separate aligned reads file into repeat reads and unique reads files             #
# (6) Make BED files from unique reads files and zip the unique reads file             #
# (7) Zip Non-aligned and Repeat files                                                 #
# (8) Convert BED files to FPKM WIG files                                              #
# (9) Convert FPKM WIG files to Visualize FPKM WIG files                               #
########################################################################################


# (5) Separate repeats from uniques into different files
SeqProcess::separate_repeats($ExperimentTopDir, $FilePrefix, $alignedpreseparationfile);

# (6) Make BED files from the unique reads bowtie output and zip the unique reads file
SeqProcess::elandext_to_bed($uniqalignedreadsfile, $ExperimentTopDir, $FilePrefix, $ReadLength, $FinalReadLength, 2, 3, 1, $MaxDupReads);
print "Zipping unique reads files\n";
`gzip $uniqalignedreadsfile`;

# (7) Zip non-aligned reads and repeat reads files 
print "Zipping non-aligned reads and repeat reads files\n";
`gzip $nonalignedreadsfile`;
`gzip $repalignedreadsfile`;

# (8) Change BED file read length (Choose the final read length)
# Create folder named "$FilePrefix_bed" inside $ExperimentTopDir to contain new bed files
#$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_extendedbed\n";
#`$commandline`;
# The original bed files contain the prefix, $FilePrefix/$FilePrefix_chr
# The new bed files (in $ExperimentTopDir) contain the prefix, $FilePrefix_bed/$FilePrefix_Chr
#my $origlengthbedfiles =  $ExperimentTopDir . $FilePrefix . "_bed/" . $FilePrefix . "_chr";
#my $finallengthbedfiles = $ExperimentTopDir . $FilePrefix . "_extendedbed/" . $FilePrefix;
#my $origlengthbedfolder =  $ExperimentTopDir . $FilePrefix . "_bed/";
#SeqProcess::change_bed_read_length($origlengthbedfiles, $finallengthbedfiles, $FinalReadLength, @Chromosomes);
# (9) Remove Unextended Bed Folder
#print "Removing unextended bed files\n";
#`rm -R $origlengthbedfolder`;


# (8) Convert BED files to FPKM WIG files
# Create folder named "$FilePrefix_FPKMWIG" inside $ExperimentTopDir to contain new FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_FPKMWIG\n";
`$commandline`;
# The WIG Track Name Prefix is $FilePrefix
SeqProcess::beddir_to_fpkmwig($bedtowigfiles, $fpkmwigfiles, $FilePrefix, $WIGTrackColor,$FinalReadLength, $MaxDupReads);

# (9) Convert FPKM WIG files to Visualize FPKMWIG files
# Create folder named "$FilePrefix_VisFPKMWIG" to contain Visualize FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_VisFPKMWIG\n";
`$commandline`;
# Step Size is set to 10 and the WIG Track Name Prefix is $FilePrefix
SeqProcess::visualize_fpkmwig($fpkmwigfiles, $visfpkmwig, 10, $WIGTrackColor, $FilePrefix);


__END__
usage: RPKM.pl 
    1) Input VarStepWIG file (format: position [tab] height)
    2) Output FPKM wig file 
    3) Read Length
    4) Read Count (in millions)

# Make RPKM for Genes

# Make RPKM for CpG Islands

__END__

my $inputbedprefix = shift(@ARGV);
my $outputbedprefix = shift(@ARGV);
my $readlengthextension = shift(@ARGV);


my @Chr;             # array that contains all the the names of the mouse chromosomes
for (my $n = 1; $n< 20; $n++){
    push(@Chr, $n);
}
push(@Chr, "M");
push(@Chr, "X");
push(@Chr, "Y");



###################
# Main Input loop #
###################

while(@Chr){
	my $chr = shift(@Chr);
	my $inputfile = $inputbedprefix . $chr . ".bed";
	open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";
	my $outfile = $outputbedprefix . "_Chr" . $chr . ".bed";
	open(OUT, ">$outfile") or die "cannot open $outfile outfile";

	print "Processing $inputfile \n";

	while(<IN>){
		chomp;
    	my @line = split ("\t", $_);
		my $start = $line[1];
		my $end = $line[2];
		if($line[5] eq "+"){
			my $end = $line[2]+$readlengthextension;
			print OUT $line[0],"\t",$line[1],"\t",$end,"\t",$line[3],"\t",$line[4],"\t",$line[5],"\n";
		}
		elsif($line[5] eq "-"){
			my $start = $line[1]-$readlengthextension;
			print OUT $line[0],"\t",$start,"\t",$line[2],"\t",$line[3],"\t",$line[4],"\t",$line[5],"\n";
		}
		else {die "$line[5] does not equal + or - \n";}
	}
}


