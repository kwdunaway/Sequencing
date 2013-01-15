#!/usr/bin/perl
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
    3) Bowtie output prefix (will have 3 files with _Uniq, _Repeat, and _Nonaligned and located in Experiment Top Folder path)
    4) Bed file prefix
    5) Read length
    6) Final read length
    7) WIG Track Color (in RRR,GGG,BBB format)
    8) Read Count (typically 1)
" unless @ARGV == 8;

my $ExperimentTopDir = shift(@ARGV);
my $rawfqfolder = shift(@ARGV);
my $BowtiePrefix = shift(@ARGV);
my $BedFilePrefix = shift(@ARGV);
my $ReadLength = shift(@ARGV);
my $FinalReadLength = shift(@ARGV);
my $WIGTrackColor = shift(@ARGV);
my $ReadCount = shift(@ARGV);

my $commandline = ""; #inputs for command line

#TODO: Must figure out how to get read count automatically
my $addtoPATH = "/home/kwdunaway/tuxedo/bowtie-0.12.7/:/home/kwdunaway/tuxedo/tophat-1.4.1.Linux_x86_64/:/home/kwdunaway/tuxedo/samtools/:/home/kwdunaway/tuxedo/cufflinks-1.3.0.Linux_x86_64/";
my $mm9path = "/home/kwdunaway/mm9/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome";


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
my $filtered_fastq = SeqProcess::filter_zip($rawfqfolder);

# (3) Run Bowtie (Separates reads between aligned and non-aligned)
my ($nonalignedreadsfile, $alignedpreseparationfile) = SeqProcess::run_bowtie($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq);

# (4) Remove made files
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
my ($uniqalignedreadsfile, $repalignedreadsfile) = SeqProcess::separate_repeats($ExperimentTopDir, $BowtiePrefix, $alignedpreseparationfile);

# (6) Zip non-aligned reads and repeat reads files 
`gzip $nonalignedreadsfile`;
`gzip $repalignedreadsfile`;


# (7) Make BED files from the Unique reads bowtie output and zip the Unique reads file
SeqProcess::elandext_to_bed($uniqalignedreadsfile, $BedFilePrefix, $ReadLength, 2, 3, 1);
`gzip $uniqalignedreadsfile`;

# (8) Change BED file read length (Choose the final read length)

# Create folder named "$BedFilePrefix_bed" inside $ExperimentTopDir to contain new bed files
$commandline = "mkdir " . $ExperimentTopDir . $BedFilePrefix . "_bed\n";
`$commandline`;

# The original bed files contain the prefix, $BedFilePrefix/$BedFilePrefix_chr
# The new bed files (in $ExperimentTopDir) contain the prefix, $BedFilePrefix_bed/$BedFilePrefix_Chr
my $origlengthbedfiles =                      $BedFilePrefix .     "/" . $BedFilePrefix . "_chr";
my $finallengthbedfiles = $ExperimentTopDir . $BedFilePrefix . "_bed/" . $BedFilePrefix;

SeqProcess::change_bed_read_length($origlengthbedfiles, $finallengthbedfiles, $FinalReadLength);

# (9) Remove Unextended Bed Folder
`rm -R $BedFilePrefix`;


# (10) BED to FPKM WIG files

# Create folder named "$BedFilePrefix_FPKMWIG" inside $ExperimentTopDir to contain new FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $BedFilePrefix . "_FPKMWIG\n";
`$commandline`;

# The bed files contain the prefix, $BedFilePrefix_bed/$BedFilePrefix_Chr
# The new FPKM WIG files contain the prefix, $BedFilePrefix_FPKMWIG/$BedFilePrefix_FPKM
my $bedtowigfiles = $ExperimentTopDir . $BedFilePrefix .     "_bed/" . $BedFilePrefix . "_Chr";
my $fpkmwigfiles =  $ExperimentTopDir . $BedFilePrefix . "_FPKMWIG/" . $BedFilePrefix . "_FPKM";

SeqProcess::vswig_to_fpkmwig($bedtowigfiles, $fpkmwigfiles, $BedFilePrefix, $WIGTrackColor, 		$FinalReadLength, $ReadCount);

# (11) Visualize FPKMWIG

# Create folder named "$BedFilePrefix_VisFPKMWIG" to contain Visualize FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $BedFilePrefix . "_VisFPKMWIG\n";
`$commandline`;

# The FPKM WIG files contain the prefix, $BedFilePrefix_FPKMWIG/$BedFilePrefix_FPKM_Chr
# The Visualize FPKM WIG files contain the prefix, $BedFilePrefix_VisFPKMWIG/$BedFilePrefix_VisFPKMWIG/
my $pre_visfpkmwig = $ExperimentTopDir . $BedFilePrefix .    "_FPKMWIG/" . $BedFilePrefix . 				"_FPKM" . "_Chr";
my $visfpkmwig =     $ExperimentTopDir . $BedFilePrefix . "_VisFPKMWIG/" . $BedFilePrefix . 				"_VisFPKMWIG/";

SeqProcess::visualize_fpkmwig($pre_visfpkmwig, $visfpkmwig, 10, $WIGTrackColor, $BedFilePrefix);


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


