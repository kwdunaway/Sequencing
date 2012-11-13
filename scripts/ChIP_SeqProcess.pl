#!/usr/bin/perl
use strict; use warnings;
use SeqProcess;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 10-30-2012
# Script Name: ChIP_SeqProcess.pl
#
# Takes a folder of bed files and creates a new one with read length extended directionally
# based on arguments input on command line.
#
# Arguments:
#    1) Input Bed File prefix (ex: DY_Chr)
#    2) Output Bed prefix without Chr (ex: NewDY will make NewDY_Chr*.bed files) 
#    3) Read Length extension (ex: 97)
#
################################################################################################

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "ChIPseq_pipeline.pl needs the following parameters:
    1) Output files name (will add _bowtie.bash and _postbowtie.bash)
    2) Experiment Top Folder path
    3) Raw file folder (make sure they are the only zipped files and the extension is .fq.gz) 
    4) Bowtie output prefix (will have 3 files with _Uniq, _Repeat, and _Nonaligned and located in Experiment Top Folder path)
    5) Bed file prefix
    6) Read length
    7) Extended read length
    8) WIG Track Color (in RRR,GGG,BBB format)
" unless @ARGV == 8;

my $outprefix = shift(@ARGV);
my $bowtieout = $outprefix . "_bowtie.bash";
open(BOWTIE, ">$bowtieout") or die "cannot open $bowtieout outfile";
my $postbowtieout = $outprefix . "_postbowtie.bash";
open(POSTBOWTIE, ">$postbowtieout") or die "cannot open $postbowtieout outfile";
my $ExperimentTopDir = shift(@ARGV);
my $rawfqfolder = shift(@ARGV);
my $BowtiePrefix = shift(@ARGV);
my $BedFilePrefix = shift(@ARGV);
my $ReadLength = shift(@ARGV);
my $ExtendedReadLength = shift(@ARGV);
my $WIGTrackColor = shift(@ARGV);

my $commandline; #inputs for command line

#TODO: Must figure out how to get read count automatically
my $ReadCount = 1;
my $addtoPATH = "/home/kwdunaway/tuxedo/bowtie-0.12.7/:/home/kwdunaway/tuxedo/tophat-1.4.1.Linux_x86_64/:/home/kwdunaway/tuxedo/samtools/:/home/kwdunaway/tuxedo/cufflinks-1.3.0.Linux_x86_64/";
my $mm9path = "/home/kwdunaway/mm9/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome";


#########################
# Bowtie and Pre-bowtie #
#########################


# Header and path modifier
SeqProcess::add_path($addtoPATH);

# Unzip zipped files, filter them, and combine into one .fq file
my $filtered_fastq = SeqProcess::filter_zip($rawfqfolder);

# Run Bowtie
my ($nonalignedreadsfile, $alignedpreseparationfile) = SeqProcess::run_bowtie($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq);

# Remove made files
`rm $filtered_fastq`;


###############
# Post-bowtie #
###############


# Separate Repeats from Uniqs
my ($uniqalignedreadsfile, $repalignedreadsfile) = SeqProcess::separate_repeats($ExperimentTopDir, $BowtiePrefix, $alignedpreseparationfile);

# Zip Nonaligned and Repeat files 
`gzip $nonalignedreadsfile`;
`gzip $repalignedreadsfile`;


# Make BED files from Uniq bowtie output and zip the Uniq file
SeqProcess::elandext_to_bed($uniqalignedreadsfile, $BedFilePrefix, $ReadLength, 2, 3, 1);
`gzip $uniqalignedreadsfile`;

# Extend BED file read length
$commandline = "mkdir "  . $ExperimentTopDir . $BedFilePrefix . "_bed\n";
`$commandline`;
my $unextendedbedfiles = $BedFilePrefix . "/" . $BedFilePrefix . "_chr";
my $extendedbedfiles = $ExperimentTopDir . $BedFilePrefix . "_bed/" . $BedFilePrefix;

SeqProcess::extend_bed_read_length($unextendedbedfiles, $extendedbedfiles, $ExtendedReadLength);

# Remove Unextended Bed Folder
`rm -R $BedFilePrefix`;


# Make WIG files
$commandline = "mkdir "  . $ExperimentTopDir . $BedFilePrefix . "_VarStepWIG\n";
`$commandline`;
print POSTBOWTIE "perl /home/kwdunaway/perl_script/BedDir_to_VarStepWIG.pl "  , $ExperimentTopDir , $BedFilePrefix , "_bed/" , $BedFilePrefix , "_Chr ", $ExperimentTopDir , $BedFilePrefix , "_VarStepWIG/", $BedFilePrefix , " " , $BedFilePrefix , " " , $WIGTrackColor , "\n";

print POSTBOWTIE "mkdir " , $ExperimentTopDir , $BedFilePrefix , "_FPKMWIG\n"; 
my $FullReadLength = $ReadLength + $ExtendedReadLength;
print POSTBOWTIE "perl /home/kwdunaway/perl_script/VarStepWIG_to_FPKMWIG.pl " , $ExperimentTopDir , $BedFilePrefix , "_VarStepWIG/", $BedFilePrefix ,  "_Chr " , $ExperimentTopDir , $BedFilePrefix , "_FPKMWIG/" , $BedFilePrefix , "_FPKM " , $FullReadLength , " " , $ReadCount , "\n";


__END__
useage: RPKM.pl 
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


