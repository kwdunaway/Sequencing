#!/usr/bin/perl
BEGIN {push @INC, "/home/kwdunaway/perl_script";}
use strict; use warnings;
use SeqProcesstemp;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 2-12-2013
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
my $test_presep = shift(@ARGV);
my $FilePrefix = shift(@ARGV);
my $FinalReadLength = shift(@ARGV);
my $WIGTrackColor = shift(@ARGV);
my $MaxDupReads = shift(@ARGV);

# Error Checking
die "Error: 6) Maximum duplicate reads should be at least 1.\n" unless $MaxDupReads > 0;

# Adjusting Folder Names
$ExperimentTopDir = $ExperimentTopDir . "\/" unless $ExperimentTopDir =~ m/\/$/;
my $commandline = ""; #inputs for command line


###########################
# Computer specific paths #
###########################

my $addtoPATH = "/home/kwdunaway/tuxedo/bowtie-0.12.7/:/home/kwdunaway/tuxedo/tophat-1.4.1.Linux_x86_64/:/home/kwdunaway/tuxedo/samtools/:/home/kwdunaway/tuxedo/cufflinks-1.3.0.Linux_x86_64/";
my $mm9path = "/home/kwdunaway/mm9/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome";
my $hg18path = "/home/kwdunaway/hg18/hg18";


####################
# Global Variables #
####################

my $ReadLength = 50;
my $nonalignedreadsfile = $ExperimentTopDir . $FilePrefix . "_NonAligned.fq";
my $alignedpreseparationfile = $test_presep;
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
SeqProcesstemp::add_path($addtoPATH);

########################################################################################
#                                Post-bowtie                                           #
#                                                                                      #
# (5) Separate aligned reads file into repeat reads and unique reads files             #
# (6) Zip Non-aligned and Repeat files                                                 #
# (7) Make BED files from unique reads files and zip the unique reads file             #
# (10) Convert BED files to FPKM WIG files                                             #
# (11) Convert FPKM WIG files to Visualize FPKM WIG files                              #
########################################################################################


# (5) Separate repeats from uniques into different files
SeqProcesstemp::separate_repeats($ExperimentTopDir, $FilePrefix, $alignedpreseparationfile);

# (7) Make BED files from the unique reads bowtie output and zip the unique reads file
SeqProcesstemp::elandext_to_bed($uniqalignedreadsfile, $ExperimentTopDir, $FilePrefix, $ReadLength, $FinalReadLength, 2, 3, 1, $MaxDupReads);
print "Zipping unique reads files\n";
`gzip $uniqalignedreadsfile`;

# (6) Zip non-aligned reads and repeat reads files 
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


# (10) BED to FPKM WIG files
# Create folder named "$FilePrefix_FPKMWIG" inside $ExperimentTopDir to contain new FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_FPKMWIG\n";
`$commandline`;
# The bed files contain the prefix, $FilePrefix_bed/$FilePrefix_Chr
# The new FPKM WIG files contain the prefix, $FilePrefix_FPKMWIG/$FilePrefix_FPKM
SeqProcesstemp::beddir_to_fpkmwig($bedtowigfiles, $fpkmwigfiles, $FilePrefix, $WIGTrackColor,$FinalReadLength, $MaxDupReads);

# (11) Visualize FPKMWIG
# Create folder named "$FilePrefix_VisFPKMWIG" to contain Visualize FPKM WIG files
$commandline = "mkdir " . $ExperimentTopDir . $FilePrefix . "_VisFPKMWIG\n";
`$commandline`;
# The FPKM WIG files contain the prefix, $FilePrefix_FPKMWIG/$FilePrefix_FPKM_Chr
# The Visualize FPKMWIG files contain the prefix, $FilePrefix_VisFPKMWIG/$FilePrefix_VisFPKMWIG_Chr
SeqProcesstemp::visualize_fpkmwig($fpkmwigfiles, $visfpkmwig, 10, $WIGTrackColor, $FilePrefix);
