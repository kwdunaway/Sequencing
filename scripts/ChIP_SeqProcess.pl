#!/usr/bin/perl
BEGIN {push @INC, "/home/kwdunaway/perl_script";}
use strict; use warnings;
use SeqProcess;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 10-30-2012
# Script Name: ChIP_SeqProcess.pl
#
# This script processes the subroutines in SeqProcess.pm.
#
# Arguments: See Below
#
################################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "ChIP_SeqProcess.pl needs the following parameters:
    1) Experiment Top Folder path (e.g. Folder/)
    2) Raw file folder (make sure they are the only zipped files and the extension is .fq.gz) 
    3) File prefix (will have 3 files with _Uniq, _Repeat, and _Nonaligned and located in Experiment Top Folder path)
    4) Final Read length (to be extended to)
    5) WIG Track Color (in RRR,GGG,BBB format)
    6) Maximum Duplicate Reads (1 for no duplicates)
" unless @ARGV == 6;

my $ExperimentTopDir = shift(@ARGV);
my $rawfqfolder = shift(@ARGV);
my $FilePrefix = shift(@ARGV);
my $FinalReadLength = shift(@ARGV);
my $WIGTrackColor = shift(@ARGV);
my $MaxDupReads = shift(@ARGV);

#Error Checking
die "Error: 1) Experiment Top Folder name should contain a '/' at the end.\n" unless $ExperimentTopDir =~ m/\/$/;
die "Error: 2) Raw File Folder name should contain a '/' at the end.\n" unless $rawfqfolder =~ m/\/$/;
die "Error: 6) Maximum duplicate reads should be at least 1.\n" unless $MaxDupReads > 0;

my $commandline = ""; #inputs for command line

#TODO: Must figure out how to get read count automatically
my $addtoPATH = "/home/kwdunaway/tuxedo/bowtie-0.12.7/:/home/kwdunaway/tuxedo/tophat-1.4.1.Linux_x86_64/:/home/kwdunaway/tuxedo/samtools/:/home/kwdunaway/tuxedo/cufflinks-1.3.0.Linux_x86_64/";
my $mm9path = "/home/kwdunaway/mm9/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome";


###########################################################################
#                          Bowtie and Pre-bowtie                          #
#                                                                         #
# (1) Pathing                                                             #
# (2) Prepare temp files for bowtie                                       #
#   (2.5) Determine Read Length                                           #
# (3) Bowtie (Separates reads between aligned and non-aligned)            #
# (4) Remove temp files                                                   #
###########################################################################

# (1) Header and path modifier
SeqProcess::add_path($addtoPATH);

# (2) Unzip zipped files, filter them, and combine into one .fq file
my $filtered_fastq = SeqProcess::filter_zip($rawfqfolder);

#   (2.5) Determine Read Length from Fastq File
my $ReadLength = SeqProcess::fastq_readlength($filtered_fastq);

# (3) Make folder for experiment and run Bowtie (separates reads between aligned and non-aligned)
my ($nonalignedreadsfile, $alignedpreseparationfile) = SeqProcess::run_bowtie($ExperimentTopDir, $FilePrefix, $mm9path, $filtered_fastq);

# (4) Remove made files
print "Removing $filtered_fastq\n";
`rm $filtered_fastq`;


########################################################################################
#                                Post-bowtie                                           #
#                                                                                      #
# (5) Separate aligned reads file into repeat reads and unique reads files             #
# (6) Zip Non-aligned and Repeat files                                                 #
# (7) Make BED files from unique reads files and zip the unique reads file             #
# (8) Extend BED file read length                                                      #
# (9) Remove unextended bed files                                                      #
# (10) Convert BED files to FPKM WIG files                                             #
# (11) Convert FPKM WIG files to Visualize FPKM WIG files                              #
########################################################################################


# (5) Separate repeats from uniques into different files
my ($uniqalignedreadsfile, $repalignedreadsfile) = SeqProcess::separate_repeats($ExperimentTopDir, $FilePrefix, $alignedpreseparationfile);

# (6) Zip non-aligned reads and repeat reads files 
print "Zipping non-aligned reads and repeat reads files\n";
`gzip $nonalignedreadsfile`;
`gzip $repalignedreadsfile`;


# (7) Make BED files from the unique reads bowtie output and zip the unique reads file
SeqProcess::elandext_to_bed($uniqalignedreadsfile, $ExperimentTopDir, $FilePrefix, $ReadLength, 2, 3, 1, $MaxDupReads);
print "Zipping unique reads files\n";
`gzip $uniqalignedreadsfile`;

# (8) Change BED file read length (Choose the final read length)

# Create folder named "$FilePrefix_bed" inside $ExperimentTopDir to contain new bed files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_extendedbed\n";
`$commandline`;

# The original bed files contain the prefix, $FilePrefix/$FilePrefix_chr
# The new bed files (in $ExperimentTopDir) contain the prefix, $FilePrefix_bed/$FilePrefix_Chr
my $origlengthbedfiles =  $ExperimentTopDir . $FilePrefix . "_bed/" . $FilePrefix . "_chr";
my $finallengthbedfiles = $ExperimentTopDir . $FilePrefix . "_extendedbed/" . $FilePrefix;

SeqProcess::change_bed_read_length($origlengthbedfiles, $finallengthbedfiles, $FinalReadLength);

# (9) Remove Unextended Bed Folder
print "Removing unextended bed files\n";
`rm -R $origlengthbedfiles`;


# (10) BED to FPKM WIG files

# Create folder named "$FilePrefix_FPKMWIG" inside $ExperimentTopDir to contain new FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_FPKMWIG\n";
`$commandline`;

# The bed files contain the prefix, $FilePrefix_bed/$FilePrefix_Chr
# The new FPKM WIG files contain the prefix, $FilePrefix_FPKMWIG/$FilePrefix_FPKM
my $bedtowigfiles = $ExperimentTopDir . $FilePrefix . "_extendedbed/" . $FilePrefix . "_Chr";
my $fpkmwigfiles =  $ExperimentTopDir . $FilePrefix . "_FPKMWIG/" . $FilePrefix . "_FPKM";

SeqProcess::vswig_to_fpkmwig($bedtowigfiles, $fpkmwigfiles, $FilePrefix, $WIGTrackColor, 		$FinalReadLength, $MaxDupReads);

# (11) Visualize FPKMWIG

# Create folder named "$FilePrefix_VisFPKMWIG" to contain Visualize FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_VisFPKMWIG\n";
`$commandline`;

# The FPKM WIG files contain the prefix, $FilePrefix_FPKMWIG/$FilePrefix_FPKM_Chr
# The Visualize FPKMWIG files contain the prefix, $FilePrefix_VisFPKMWIG/$FilePrefix_VisFPKMWIG_Chr
my $pre_visfpkmwig = $ExperimentTopDir . $FilePrefix .    "_FPKMWIG/" . $FilePrefix . 				"_FPKM" . "_Chr";
my $visfpkmwig =     $ExperimentTopDir . $FilePrefix . "_VisFPKMWIG/" . $FilePrefix . 				"_VisFPKMWIG";

SeqProcess::visualize_fpkmwig($pre_visfpkmwig, $visfpkmwig, 10, $WIGTrackColor, $FilePrefix);


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


