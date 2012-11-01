package SeqProcess;
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 10-30-2012
# Module Name: SeqProcess.pm
#
# This is a module with sequencing processing commands.
#
################################################################################################

###########################################################################
#                        Header and Path Modifier                         #
###########################################################################

sub add_path {
	my ($addtoPATH) = @_;

	`#!/bin/bash\n\n`;
	`PATH=\$PATH:$addtoPATH\n`;
	`export PATH\n\n`;
}

###########################################################################
#                 Combine and Filter Zipped Files                         #
#  Input: Raw file folder (only zipped files and the extension is .fq.gz) #
# Output: Filtered and Combined into one .fq file                         #
###########################################################################

sub filter_zip {
	my ($rawfqfolder) = @_;
	my $filtered_fastq = $rawfqfolder . "filtered.fq";

	`gunzip -c $rawfqfolder *.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v \"^--\$\" >  $filtered_fastq\n\n`;

	return $filtered_fastq;
}
	
###########################################################################
#                               Run Bowtie                                #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) MM9 Path                                                     #
#         4) Filtered Fastq File                                          #
#                                                                         #
# Output: 1) Non-aligned Reads File                                       #
#         2) Aligned Preseparation File                                   #
###########################################################################

sub run_bowtie {
	my ($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq) = @_;
	
	`mkdir $ExperimentTopDir\n`;

	my $nonalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_NonAligned.fq";
	my $alignedpreseparationfile = $ExperimentTopDir . $BowtiePrefix . "_alignedpreseparation.txt";

	`bowtie -p 4 -M 1 -k 1 --chunkmbs 256 --strata --best --un $nonalignedreadsfile $mm9path $filtered_fastq $alignedpreseparationfile\n\n`;

	return ($nonalignedreadsfile, $alignedpreseparationfile);
}

###########################################################################
#                     Separate Repeats from Uniques                       #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) Aligned Preseparation File                                   #
#                                                                         #
# Output: 1) Unique Aligned Reads File                                    #
#         2) Repetitive Aligned Reads File                                #
###########################################################################

sub separate_repeats {
	my ($ExperimentTopDir, $BowtiePrefix, $alignedpreseparationfile) = @_;

	my $uniqalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Uniq.txt";
	my $repalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Repeat.txt";

	open(IN, "<$alignedpreseparationfile") or die "cannot open $alignedpreseparationfile infile";
	open(UNIQOUT, ">$uniqalignedreadsfile") or die "cannot open $uniqalignedreadsfile outfile";
	open(REPOUT, ">$repalignedreadsfile") or die "cannot open $repalignedreadsfile outfile";

	while (<IN>)
	{
		chomp;
		my @line = split ("\t", $_);
		if($line[6] > 0)
		{
			print REPOUT $_, "\n";
		}
		else
		{
			print UNIQOUT $_, "\n";
		}
	}


	close IN;
	close UNIQOUT;
	close REPOUT;

	return ($uniqalignedreadsfile, $repalignedreadsfile);
}

1;

__END__

###########################################################################
#                     Eland Extended Format to BED                        #
#                                                                         #
#  This subroutine takes sequences in Eland Extended format (also known   #
#  as s_#_export.txt files) and produces multiple BED files which can     #
#  be used to analyze the data.                                           #
#                                                                         #
#  Input: 1) Input file name                                              #
#         2) Output file name                                             #
#         3) Read length                                                  #
#         4) Chromosome array number                                      #
#         5) Position array number                                        #
#         6) Strand array number                                          #
#         7) First Char                                                   #
#                                                                         #
# Output: 1)                                                              #
###########################################################################

sub elandext_to_bed {
	my ($infile, $outfile, $readlength, $chr, $pos, $strand, $firstchar) = @_;

	open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read

	my @array;
	my $QCcount = 0;
	my $NMcount = 0;
	my $NonUniqcount = 0;
	my $Unknowncount = 0;
	my @ChrMapcount;
	for (my $n = 0; $n < 27; $n++)
	{
		$ChrMapcount[$n] = 0;
	}
	my $totalcount = 0;
	my $chrnum = 1;
	my $filename;

	for(

	$filename = $outfile . "/" . $outfile . "_chrX" . ".bed";
	open(OUTX, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
  
	$filename = $outfile . "/" . $outfile . "_chrY" . ".bed";
	open(OUTY, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

	$filename = $outfile . "/" . $outfile . "_chrM" . ".bed";
	open(OUTM, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written
  
	$filename = $outfile . "/" . $outfile . "_NM.fq";
	open(OUTNM, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

	$filename = $outfile . "/" . $outfile . "_NonUnique.fq";
	open(OUTNONUNIQ, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

	$filename = $outfile . "/" . "Stats_" . $outfile . ".txt";
	open(OUTSTAT, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

	$filename = $outfile . "/" . "Unknown_" . $outfile . ".txt";
	open(UNKNOWN, ">$filename") or die "cannot open $filename outfile"; #opens outfile to be written

}
# BED to WIG
# WIG to FPKMWIG
# Vis FPKMWIG



1;

