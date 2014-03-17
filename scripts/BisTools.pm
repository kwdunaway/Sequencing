package BisTools;
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 12-30-2012
# Module Name: BisTools.pm
#
# This is a module with sequencing processing commands related to aligning and
# SPKM (Segments Per Kilobase per Million fragments mapped).
#
# Jump to each subroutine by searching for the number (e.g. "(6)" ) or name (e.g. "Sort BED File" ).
# Descriptions of each subroutine can be found above their code.
#
# Subroutines: 
#              (0) Determine Read Length from Fastq.gz File
#                (0.1) Determine Read Length from Fastq File
#                (0.2) Determine Read Count from Fastq File
#                (0.3) Determine Read Count from Eland Extended File
#              (1) Header and Path Modifier
#
#              -Bowtie Related-
#              (2) Combine and Filter Zipped Fastq Files
#              (3) Run Bowtie for Non-aligned and Aligned Reads
#              (4) Separate Aligned Reads - Repeats from Uniques
#
#              -Eland Extended/BED Related-
#              (5) Eland Extended Format to BED
#              (6) Sort BED File
#              (7) Eliminate Duplicate Reads in BED Files
#              (8) Extend Read Length of BED Files
#              (9) Change Read Length of BED Files (Choose Read Length)
#
#              -RPKM/FPKM/WIG Related-
#             (10) BED Directory to Variable Step WIG
#             (11) Variable Step WIG to FPKM WIG
#             (12) BED Directory to FPKM WIG
#             (13) Visualize FPKM WIG
#             (14) RPKM From BED and GTF
#             (15) FPKM from GTF and FPKMWIG
#
##########################################################################################

###########################################################################
#       (0) Determine Read Length from Gzipped Fastq (.fq.gz) File        #
###########################################################################

sub fastqgz_readlength 
{
	my ($fastqgzfile) = @_; # Input is "filename.fq.gz"
	print "Read Length from .fq.gz File: Obtaining read length from " , $fastqgzfile , "\n";
	# Uses gunzip to determine read length from the second line of the file
	my $readlength = `gunzip -c $fastqgzfile | head -n 2 | tail -n 1 | tr -d '\n'| wc -m | tr -d '\n'`;
	print "Read Length from .fq.gz File: Finished. Read length is ", $readlength ," \n\n";
	return $readlength;
}

