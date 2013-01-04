package SeqProcess;
use strict; use warnings;

###############################################################################################
# Author: Keith Dunaway & Roy Chu
# Email: kwdunaway@ucdavis.edu & rgchu@ucdavis.edu
# Date: 12-30-2012
# Module Name: SeqProcess.pm
#
# This is a module with sequencing processing commands related to aligning and
# SPKM (Segments Per Kilobase per Million fragments mapped).
#
# Jump to each subroutine by searching for the number (e.g. "(6)" ) or name (e.g. "Sort BED File" ).
# Descriptions of each subroutine can be found above their code.
#
# Subroutines: 
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
################################################################################################

###########################################################################
#                     (1) Header and Path Modifier                        #
###########################################################################

sub add_path 
{
	my ($addtoPATH) = @_;

	`#!/bin/bash\n\n`;
	`PATH=\$PATH:$addtoPATH\n`;
	`export PATH\n\n`;
}

###########################################################################
#              (2) Combine and Filter Zipped Fastq Files                  #
#  Input: Raw file folder (only zipped files and the extension is .fq.gz) #
# Output: Returns: Filtered and Combined into one .fq file                #
###########################################################################

sub filter_zip 
{
	my ($rawfqfolder) = @_;
	my $filtered_fastq = $rawfqfolder . "filtered.fq";

	`gunzip -c $rawfqfolder*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v \"^--\$\" >  $filtered_fastq\n\n`;

	return $filtered_fastq;
}
	
###########################################################################
#            (3) Run Bowtie for Non-aligned and Aligned Reads             #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) MM9 Path                                                     #
#         4) Filtered Fastq File                                          #
#                                                                         #
# Output: Returns: 1) Non-aligned Reads File                              #
#                  2) Aligned Preseparation File                          #
###########################################################################

sub run_bowtie 
{
	my ($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq) = @_;
	
	`mkdir $ExperimentTopDir\n`;

	my $nonalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_NonAligned.fq";
	my $alignedpreseparationfile = $ExperimentTopDir . $BowtiePrefix . "_alignedpreseparation.txt";

	`bowtie -p 4 -M 1 -k 1 --chunkmbs 256 --strata --best --un $nonalignedreadsfile $mm9path $filtered_fastq $alignedpreseparationfile\n\n`;

	return ($nonalignedreadsfile, $alignedpreseparationfile);
}

###########################################################################
#             (4) Separate Aligned Reads - Repeats from Uniques           #
#  Input: 1) Experiment Top Folder path                                   #
#         2) Bowtie output prefix                                         #
#         3) Aligned Preseparation File                                   #
#                                                                         #
# Output: Returns: 1) Unique Aligned Reads File                           #
#                  2) Repetitive Aligned Reads File                       #
###########################################################################

sub separate_repeats 
{
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

###########################################################################
#                   (5) Eland Extended Format to BED                      #
#   *             Note: subroutine calls the "Sort BED"               *   #
#   *        and "Eliminate Duplicate BED Reads" subroutines          *   #
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
#                                                                         #
# Output: Creates directory containing:                                   #
#           1) BED files for each chromosome                              #
#           2) Statistics on all reads                                    #
#           3) Gzipped non-unique and non-mappable reads files            #
#           4) Unknown reads file                                         #
###########################################################################

sub elandext_to_bed 
{
	# Input
	my ($infile, $outfile, $readlength, $chr, $pos, $strand) = @_;

	# Makes Output Directory
	if (! -d $outfile) 
	{ 
		`mkdir $outfile`; #creates dir if one doesn't exist
		if (! -d $outfile) { die "directory ($outfile) does not exist"} 
	}

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be read

	my @array; # Contains the data to be parsed
	my $QCcount = 0; # Low Quality Reads (placed in non-mappable reads)
	my $NMcount = 0; # Non-mappable Reads
	my $NonUniqcount = 0; # Non-unique Reads
	my $Unknowncount = 0; # Unknown Reads
	my @ChrMapcount; # Unique Mapped Reads by Chromosome
	for (my $n = 0; $n < 27; $n++)
	{
		$ChrMapcount[$n] = 0;
	}
	my $totalcount = 0; #Total Mapped Reads
	my $filename; # Outfile for each chromosome in bed
	my @chr_out; # Holds all the output files
	my $commandinput; # Holds Command Line Inputs

	# Hash made to include X, Y, and M
	my %chromosomes = (); # Contains chromosomes(keys = call number, values = chromosome number)

	for(my $n = 1; $n < 24; $n++)
	{
		$chromosomes{$n} = $n;
	}
	$chromosomes{24} = "X";
	$chromosomes{25} = "Y";
	$chromosomes{26} = "M";

	############################################################################
	#       Creates BED file for each Chromosome, Puts Output files in array   #
	#                                                                          #
	#   Chromosome 1 corresponds to $chr_out[0], Chr2 -> $chr_out[1], etc.     #
	#                    ChrX = [23], ChrY = [24], ChrMitochondria = [25]      #
	#   *Fastq Files* -> Non-mappable = [26], Non-unique = [27]                #
	#   *Text Files*  -> Stats = [28], Unknown = [29]                          #
	############################################################################

	for(my $n = 0; $n < 30; $n++)
	{
		if($n > 22)
		{
			if($n == 23)
			{
				$filename = $outfile . "/" . $outfile . "_chrX.bed";
			}
			if($n == 24)
			{
				$filename = $outfile . "/" . $outfile . "_chrY.bed";
			}
			if($n == 25)
			{
				$filename = $outfile . "/" . $outfile . "_chrM.bed";
			}
			if($n == 26)
			{
				$filename = $outfile . "/" . $outfile . "_NM.fq"; 
			}
			if($n == 27)
			{
				$filename = $outfile . "/" . $outfile . "_NonUnique.fq";
			}
			if($n == 28)
			{
				$filename = $outfile . "/" . "Stats_" . $outfile . ".txt";
			}
			if($n == 29)
			{
				$filename = $outfile . "/" . "Unknown_" . $outfile . ".txt";
			}
		}
		else
		{
			$filename = $outfile . "/" . $outfile . "_chr" . ($n + 1) . ".bed";
		}
		open($chr_out[$n], ">$filename") or die "cannot open $filename outfile";
	}


	############################################################################
	#           Processing data into each chromosome output file               #
	############################################################################

	while (<IN>)
	{
		chomp;
		@array = split("\t", $_); # Splitting data into array
		$totalcount++; # Total mapped reads
		# Non-mappable reads
		if ($array[$chr] eq "QC") # Low quality reads (placed in non-mappable reads)
		{ 
			$chr_out[26]->print("$_" , "\n"); 
			$QCcount++;
		}
		elsif ($array[$chr] eq "NM") # Non-mappable reads
		{ 	
			$chr_out[26]->print("$_" , "\n");  
			$NMcount++;
		}
		# Non-unique reads
		elsif ($array[6] > 0)
		{ 	
			$chr_out[27]->print("$_" , "\n"); 
			$NonUniqcount++;
		}
		# Processing unique reads to chromosomes
		elsif($array[$chr] =~ m/chr/)
		{
			foreach my $chromosome(keys %chromosomes) # Check all chromosomes
			{	
				my $chr_value = $chromosomes{$chromosome};
				if($array[$chr] =~ m/(chr$chr_value)$/)
				{ 	
					$chr_out[$chromosome-1]->print("chr", $chr_value, "\t" , 
								$array[$pos] , "\t" , 
								$array[$pos]+$readlength, 									"\t", $outfile , "\t", "0", "\t" , 
								$array[$strand] , "\n");
					$ChrMapcount[$chromosome]++;
				}
			}
			$ChrMapcount[0]++; # Add to total number of uniquely mapped reads
		}
		else # Unknown data
		{
			$chr_out[29]->print("$_\n");
			$Unknowncount = $Unknowncount + 1;
  		}
	}

	############################################################################
	#                  Printing statistics to Stats Outfile                    #
	############################################################################

	my $weirdcount = $totalcount - $QCcount - $NMcount - $NonUniqcount - $ChrMapcount[0];

	for (my $n = 1; $n < 24; $n++)
	{
		if ($ChrMapcount[$n] > 0)
		{
			$chr_out[28]->print("Number of reads mapped to Chromosome $n is:\t" , 
								$ChrMapcount[$n], "\n");
		}
	}
	$chr_out[28]->print("Number of reads mapped to Chromosome X is:\t", $ChrMapcount[24], "\n");
	$chr_out[28]->print("Number of reads mapped to Chromosome Y is:\t", $ChrMapcount[25], "\n");
	$chr_out[28]->print("Number of reads mapped to Mitochondria is:\t", $ChrMapcount[26], "\n\n");
	$chr_out[28]->print("Number of Uniquely mapped reads is:\t", $ChrMapcount[0], "\n");
	$chr_out[28]->print("Number of NonUniquely mapped reads is:\t", $NonUniqcount, "\n");
	$chr_out[28]->print("Number of QC (Low Quality reads) is:\t", $QCcount, "\n");
	$chr_out[28]->print("Number of NM (nonmappable reads) is:\t", $NMcount, "\n");
	$chr_out[28]->print("Number of unknown results (printed in separate file) is:\t", 
								$Unknowncount, "\n");
	$chr_out[28]->print("Number of weird reads (should be same as previous line) is:\t", 
								$weirdcount, "\n");
	$chr_out[28]->print("Number of Total mapped reads is:\t", $totalcount, "\n");

	############################################################################
	#      Closing Files and Gzip Non-mappable and Non-unique Reads Files      #
	############################################################################

	# Closing files
	for(my $n = 0; $n < 30; $n++)
	{
		close($chr_out[$n]);
	}
	close(IN);

	# Gzip Non-mappable and non-unique reads outfiles
	$commandinput = "gzip " . $outfile . "/" . $outfile . "_NM.fq";
	`$commandinput`;
	$commandinput = "gzip " . $outfile . "/" . $outfile . "_NonUnique.fq";
	`$commandinput`;

	############################################################################
	#               Sort All Bed Files and Delete Duplicate Reads              #
        #  *Uses "Sort Bed File" and "Eliminate Duplicate BED Reads" subroutines*  #
	############################################################################

	my $bedfile;

	foreach my $chromosome(keys %chromosomes) # Check all chromosomes
	{
		my $chr_value = $chromosomes{$chromosome};
		$bedfile = $outfile . "/" . $outfile . "_chr" . $chr_value . ".bed";
		if ($ChrMapcount[$chromosome] == 0) # Check if empty
		{
			`rm $bedfile`; # Delete if empty

		}
		else
		{
			sort_bed($bedfile);
			eliminate_bed_dups($bedfile);
		}
	}
}

###########################################################################
#                           (6) Sort BED File                             #
#  Input: Unsorted BED file                                               #
# Output: Sorted BED file (replaces unsorted file)                        #
###########################################################################

sub sort_bed
{
	my ($bedfile) = @_;
	my $temp = $bedfile . "_sorted";
	`sort -n +1 -2 $bedfile > $temp`;
	`rm $bedfile`;
	`mv $temp $bedfile`;
}

###########################################################################
#               (7) Eliminate Duplicate Reads in BED Files                #
# Checks BED files for duplicate reads and deletes the duplicates         #
#                                                                         #
#  Input: Sorted BED File                                                 #
# Output: BED File cleaned of duplicates (replaces input file)            #
###########################################################################

sub eliminate_bed_dups
{
	# Input
	my ($bedfile) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################
	
	my $temp = $bedfile . "_nodups"; # Temporary output file

	open(IN, "<$bedfile") or die "cannot open $bedfile bed file";
	open(OUT, ">$temp") or die "cannot open $temp temp file";

	my $linenum = 1; # Record line number(according to output file)
	# Check if lines are the same; if so, duplicates
	my %data; # Will hold data of every line for checking

	##################################################
	#            Checking for Duplicates             #
	##################################################

	# Retrieve data for first line
	$_ = <IN>;
	$data{$linenum} = $_;
	my @line = split("\t", $_);
	print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 			$line[4], "\t", $line[5];
	$linenum++;

	while (<IN>)
	{
		$data{$linenum} = $_;
		my @line = split("\t", $_);

		# If current line does not equal previous line, print to output
		if ($data{$linenum-1} ne $data{$linenum})
		{
			print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 				$line[4], "\t", $line[5];

			$linenum++; 
		}
	}

	close(IN);
	close(OUT);

	`rm $bedfile`;
	`mv $temp $bedfile`;
}

########################################################################################
#                   (8) Extend Read Length of BED Files                                #
# Takes a folder of bed files and creates a new one with read length                   #
# extended directionally based on arguments input on command line.                     #
#                                                                                      #
#  Input: 1) Input Bed File prefix (ex: DY_Chr)                                        #
#         2) Output Bed prefix without Chr (ex: NewDY will make NewDY_Chr*.bed files)  #
#         3) Read Length extension (ex: 97)(added to old read length)                  #
#                                                                                      #
# Output: Extended Read Length BED Files                                               #
########################################################################################

sub extend_bed_read_length
{
	# Input
	my ($inputbedprefix, $outputbedprefix, $readlengthextension) = @_;

	my @Chr; 	# array that contains all the the names of the mouse chromosomes
	for (my $n = 1; $n< 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	while(@Chr)
	{
		my $chr = shift(@Chr);
		my $inputfile = $inputbedprefix . $chr . ".bed";
		open(IN, "<$inputfile") or die "cannot open $inputfile infile";
		my $outfile = $outputbedprefix . "_Chr" . $chr . ".bed";
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";

		print "Processing $inputfile \n";

		while(<IN>)
		{
			chomp;
    			my @line = split ("\t", $_);
			my $start = $line[1];
			my $end = $line[2];
			if($line[5] eq "+")
			{
				my $end = $line[2]+$readlengthextension;
				print OUT $line[0],"\t",$line[1],"\t",$end,"\t",$line[3],"\t",
					$line[4],"\t",$line[5],"\n";
			}
			elsif($line[5] eq "-")
			{
				my $start = $line[1]-$readlengthextension;
				print OUT $line[0],"\t",$start,"\t",$line[2],"\t",$line[3],"\t",
					$line[4],"\t",$line[5],"\n";
			}
			else {die "$line[5] does not equal + or - \n";}
		}

		close(IN);
		close(OUT);
	}
}

########################################################################################
#            (9) Change Read Length of BED Files (Choose Read Length)                  #
# Takes a folder of bed files and creates a new one with read length based on          #
# arguments input on command line (This is a version of extend_bed_read_length).       #
#                                                                                      #
#  Input: 1) Input Bed File prefix (ex: DY_Chr)                                        #
#         2) Output Bed prefix without Chr (ex: NewDY will make NewDY_Chr*.bed files)  #
#         3) New Read Length (ex: 97)                                                  #
#                                                                                      #
# Output: New Read Length BED Files                                                    #
########################################################################################

sub change_bed_read_length
{
	# Input
	my ($inputbedprefix, $outputbedprefix, $readlengthextension) = @_;

	my @Chr;	# array that contains all the the names of the mouse chromosomes
	for (my $n = 1; $n< 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	while(@Chr)
	{
		my $chr = shift(@Chr);
		my $inputfile = $inputbedprefix . $chr . ".bed";
		open(IN, "<$inputfile") or die "cannot open $inputfile infile";
		my $outfile = $outputbedprefix . "_Chr" . $chr . ".bed";
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";

		print "Processing $inputfile \n";

		while(<IN>)
		{
			chomp;
    			my @line = split ("\t", $_);
			my $start = $line[1];
			my $end = $line[2];
			if($line[5] eq "+")
			{
				my $end = $line[1]+$readlengthextension;
				print OUT $line[0],"\t",$line[1],"\t",$end,"\t",$line[3],"\t",
					$line[4],"\t",$line[5],"\n";
			}
			elsif($line[5] eq "-")
			{
				my $start = $line[2]-$readlengthextension;
				print OUT $line[0],"\t",$start,"\t",$line[2],"\t",$line[3],"\t",
					$line[4],"\t",$line[5],"\n";
			}
			else {die "$line[5] does not equal + or - \n";}
		}

		close(IN);
		close(OUT);
	}
}

###########################################################################
#              (10) BED Directory to Variable Step WIG                    #
# Converts a directory of BED files to a directory of WIG files           #
#                                                                         #
#  Input: 1) Input BED File prefix (ex: DY_Chr)                           #
#         2) Output WIG prefix without Chr (ex: NewDY -> NewDY_Chr*.wig)  #
#         3) WIG Track Name Prefix (for genome browser)                   #
#         4) WIG Track Color (Format: RRR,GGG,BBB)                        #
#                                                                         #
# Output: Directory of .wig Files                                         #
###########################################################################

sub beddir_to_vswig
{
	# Input
	my ($infileroot, $outfileroot, $wignameroot, $color) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	my @Chr;	# array that contains all the the names of the chromosomes

	for (my $n = 1; $n < 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	####################################################################################
	# Grabs the information from each line and assigns it to the appropriate variables #
	####################################################################################

	print "\n\nStarting Bed to Wig conversion of files with prefix $infileroot:\n";
	while(@Chr)
	{
		my $infile = $infileroot . $Chr[0] . ".bed";
		open(IN, "<$infile") or die "cannot open $infile infile"; #opens input file to be 											read (must be .bed)
		my @infileroot = split(".bed", $infile);
		my $outfile = $outfileroot . "_Chr" . $Chr[0] . ".wig";
		open(OUT, ">$outfile") or die "cannot open $outfile outfile"; #opens output file to 											write to (.wig)
		# Prints the head of the track (necessary for genome browser to read file properly) 			(customizable through terminal)
		my $wigname = $wignameroot . "_Chr" . $Chr[0];
		print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", 
			$wigname, "\" description=\"", $wigname, "\" color=", $color, "\n";

		my %PosVal;              # hash that contains the peak heights, kept small
		my $position = 0;        # position of window
		my $chromcheck;          # variable that checks each chromosome to make sure they 							are the same chromosome
		my @line;                # temp array used to retreive the information of each line
		my $chrom;               # chromosome of current line
		my $startread = 0;       # start of read at current line
		my $endread = 0;         # end of read at current line 						(does not include this position)

		while (<IN>)
		{
			chomp;
			@line = split ("\t", $_);
			$chrom = $line[0];	# makes $chrom have the chromosome information
			$startread = $line[1];	# sets the value for the start of the read to 							$startread
			$endread = $line[2];	# sets the value for the end of the read to $endread

	    		# Makes the program run faster by skipping the first gap
			if ($startread > -1)  # gets rid of all lines before 0 point
			{
				#    print "$chrom \t $startread \t $endread \n";
				if ($position == 0)
		    		{ 
					$position = $startread;	
					$chromcheck = $chrom;
					print OUT "variableStep chrom=", $chrom," step=1\n";
				}
	    
				# ends program if you have different chromosomes in your data
				die "You have different chromosomes in this bed file (Program ended before completion)/n" unless $chrom == $chromcheck;
	    
				# adds height to PosVal for the sequence found
				my $addcount = $startread;
				while ($addcount < $endread)
				{
					if (exists $PosVal{$addcount})
					{
						$PosVal{$addcount} = $PosVal{$addcount} + 1;
						# print "$PosVal{$addcount} \n";
					}
					else
					{
						$PosVal{$addcount} = 1;
					}
					# print "$addcount \t",$PosVal{$addcount}, "\n";
					$addcount = $addcount +1;
				}

	##################################################
	# Prints data to a fixed step wiggle file (.wig) #
	##################################################

	# Since you should NEVER have a read's start before the current read's start, 
	# this will print all positions until that point

				while ($position < $startread) 
				{
					if (exists $PosVal{$position})
					{
						print OUT $position, "\t", $PosVal{$position}, "\n";
						delete $PosVal{$position};
						$position = $position + 1;
					}
					else # skips gaps by jumping to the next position and 							starting a new fixed step line
					{
		    				$position = $startread;
					}
				}
			}
		}

	##############################################################
	# Prints data past the window to $outfile (ONE LAST TIME!!!) #
	##############################################################

	#  This will print the rest of the reads
		while ($position < $endread) 
		{
			if (exists $PosVal{$position})
			{
				print OUT $position, "\t", $PosVal{$position}, "\n";
				delete $PosVal{$position};
				$position = $position + 1;
			}
		}
    
		close IN;
		close OUT;
		print "Finished with Chromosome $Chr[0] \n";
		shift(@Chr);
	}

}

###########################################################################
#                  (11) Variable Step WIG to FPKM WIG                     #
# Calculates FPKM (Fragments Per Kilobase per Million mapped reads)       #
# from raw WIG files using read length and read count (in millions)       #
#                                                                         #
#  Input: 1) Input VarStepWIG prefix (format: position [tab] height)      #
#         2) Output FPKM WIG prefix without Chr                           #
#         3) Read Length                                                  #
#         4) Read Count (in millions)                                     #
#                                                                         #
# Output: Directory of FPKM WIG Files                                     #
###########################################################################

sub vswig_to_fpkmwig
{
	# Input
	my ($inputWIGprefix, $outprefix, $readlength, $readcount) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	my @Chr;             # array that contains all the the names of the mouse chromosomes
	for (my $n = 1; $n < 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	#############################################
	# Calculate and print OUT FPKM from file(s) #
	#############################################

	while(@Chr)
	{
		my $chr = shift(@Chr);
		print "Now Converting: Chr$chr\n";
		my $inputWIG = $inputWIGprefix . $chr . ".wig";
		my $outfile = $outprefix . "_Chr" . $chr . ".wig";
		open(INWIG, "<$inputWIG") or die "cannot open $inputWIG INWIG infile";
		open(OUT, ">$outfile") or die "cannot open $outfile outfile";
		my $lin = <INWIG>;
		print OUT $lin;
		$lin = <INWIG>;
		print OUT $lin;
		while(<INWIG>)
		{
			chomp;
			my @line = split ("\t", $_);
			if(exists $line[1])
			{
	    			my $FPKM = ($line[1] * 1000) / ($readlength * $readcount);
				$FPKM = sprintf("%.4f",$FPKM);
				print OUT $line[0], "\t" , $FPKM,"\n";
			}
		}
		close OUT;
		close INWIG;
	}
	
}

###########################################################################
#                    (12) BED Directory to FPKM WIG                       #
#   * Note: subroutine calls the "BED Directory to Variable Step WIG" *   #
#   *        and "Variable Step WIG to FPKM WIG" subroutines          *   #
#                                                                         #
# Converts a directory of BED files to a directory of FPKM WIG files      #
#                                                                         #
#  Input: 1) Input BED File prefix (ex: DY_Chr)                           #
#         2) Output FPKM WIG prefix without Chr                           #
#         3) WIG Track Name Prefix (for genome browser)                   #
#         4) WIG Track Color (Format: RRR,GGG,BBB)                        #
#         5) Read Length                                                  #
#         6) Read Count (in millions)                                     #
#                                                                         #
# Output: Directory of FPKM WIG Files                                     #
###########################################################################

sub beddir_to_fpkmwig
{
	# Input
	my ($inputprefix, $outprefix, $wignameroot, $color, $readlength, $readcount) = @_;

	my $varstepprefix = $inputprefix . "_TempVarStep";
	my $varstepfiles  = $inputprefix . "_TempVarStep" . "_Chr";

	# BED Directory to Variable Step WIG
	beddir_to_vswig($inputprefix, $varstepprefix, $wignameroot, $color);

	# Variable Step WIG to FPKM WIG
	vswig_to_fpkmwig($varstepfiles, $outprefix, $readlength, $readcount);

	# Remove Temporary Files
	`rm $varstepfiles*`;
}

###########################################################################
#                         (13) Visualize FPKM WIG                         #
# Combines FPKM (Fragments Per Kilobase per Million mapped reads)         #
# files given multipliers for them                                        #
#                                                                         #
#  Input: 1) Input prefix for FPKM WIG files(leave out chr# and .wig)     #
#         2) Output prefix FPKM WIG files without Chr                     #
#         3) Step Size                                                    #
#         4) WIG Track Color (Format: RRR,GGG,BBB)                        #
#         5) WIG Track Name Prefix (for genome browser)                   #
#                                                                         #
# Output: Gzipped FPKM WIG Files                                          #
###########################################################################

sub visualize_fpkmwig
{
	# Input
	my ($inputprefix, $outprefix, $stepsize, $color, $tracknameprefix) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################
	my @Chr;             # array that contains all the the names of the mouse chromosomes
	for (my $n = 1; $n < 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	#############################################
	# Visualize and print OUT FPKM from file(s) #
	#############################################

	while(@Chr)
	{
		my $chr = shift(@Chr);
		print "Now Visualizing: Chr$chr\n";
		my $inputWIG = $inputprefix . $chr . ".wig";
		open(IN, "<$inputWIG") or die "cannot open $inputWIG IN infile";
		# Shave off first two header lines
		my $lin = <IN>;
		$lin = <IN>;
		my $outfile = $outprefix . "_Chr" . $chr . ".wig";	
		open(OUT, ">$outfile") or die "cannot open $outfile OUT outfile";
		my $trackname = $tracknameprefix . "_Chr" . $chr;
	
		print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $trackname, "\" description=\"", $trackname, "\" color=", $color, "\n";
		print OUT "variableStep chrom=chr", $chr," span=", $stepsize ,"\n";

		my %HeightHash;
		$lin = <IN>;
		my @firstline = split("\t",substr($lin, 0, -1));
		my $currentpos = $firstline[0];
		my $startpos = $firstline[0];
		$HeightHash{$firstline[0]}=$firstline[1];
		while(<IN>)
		{
			chomp;
			my @line = split("\t",$_);
			$HeightHash{$line[0]}=$line[1];
			$currentpos = $line[0];
			if ($startpos + $stepsize <= $currentpos)
			{
				my $height=0;
				for(my $n = $startpos; $n < ($startpos + $stepsize); $n++)
				{
					if(exists $HeightHash{$n})
					{	
						$height = $height + $HeightHash{$n};
						delete $HeightHash{$n};
					}
				}
				$height = $height / $stepsize;
				$height = sprintf("%.2f", $height);
				print OUT $startpos , "\t" , $height , "\n";
				$startpos = $currentpos;
			}
		}
		my $height=0;
		for(my $n = $startpos; $n < ($startpos + $stepsize); $n++)
		{
			if(exists $HeightHash{$n})
			{	
				$height = $height + $HeightHash{$n};
			}
		}
		$height = $height / $stepsize;
		$height = sprintf("%.5f", $height);
		print OUT $startpos , "\t" , $height , "\n";
		close IN;
		close OUT;	
		print "Finished making $outfile, now zipping it.\n";
		my $commandline = "gzip $outfile";
		`$commandline`;
	}
}

###########################################################################
#                       (14) RPKM From BED and GTF                        #
# Scores RPKM (Reads Per Kilobase per Million mapped reads)               #
# for the given GTF file and Bed directory/prefix                         #
#                                                                         #
#  Input: 1) Input BED prefix						  #
#         2) Input GTF file                                               #
#         3) Output RPKM table file name                                  #
#         4) Total Reads (in millions)                                    #
#                                                                         #
# Output: Returns: RPKM Output File                                       #
###########################################################################

sub rpkm_from_bed
{
	# Input
	my ($inputBEDprefix, $GTFfilename, $OutputName, $TotalReads) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(GTF, "<$GTFfilename") or die "cannot open $GTFfilename GTF infile";
	open(OUT, ">$OutputName") or die "cannot open $OutputName OUT outfile";
	print OUT "Gene_Name\tChrom_pos\tGene_length\tStrand\tRPKM\n";

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

	# Takes a GTF file with lines like this (tab-delimited):
	#     chrom   source  feature  start    end    score  strand  frame  attribute
	# ex: chr1   unknown   exon   3204563  3207049   .      -       .   gene_id "Xkr4"; gene_name "Xkr4"; p_id "P2690"; transcript_id "NM_001011874"; tss_id "TSS1818";
	# and returns the following data structure:
	#
	# GTF{chrom}->{gene_name}->[start,stop,strand,count]
	# 
	# Note: count will be used later in the script but for now will be 0.

	print "Loading GTFHash from GTF file.\n";
	my %GTFHash;
	while(<GTF>)
	{
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
		if($feature eq "exon")
		{
			if(exists $GTFHash{$chrom}{$gene_name})
			{
				if($start < $GTFHash{$chrom}{$gene_name}[0])
				{
					$GTFHash{$chrom}{$gene_name}[0] = $start;
				}
				if($end < $GTFHash{$chrom}{$gene_name}[0])
				{
					$GTFHash{$chrom}{$gene_name}[0] = $end;
				}
				if($end > $GTFHash{$chrom}{$gene_name}[1])
				{
					$GTFHash{$chrom}{$gene_name}[1] = $end;
				}
				if($start > $GTFHash{$chrom}{$gene_name}[1])
				{
					$GTFHash{$chrom}{$gene_name}[1] = $start;
				}
			}
			else
			{
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

	while(@Chr)
	{
		my $chr = shift(@Chr);
		my $chrom_name = "chr" . $chr;

		#opens infile
		my $inputfile = $inputBEDprefix . $chr . ".bed";
		open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";

		#Makes a StartHash for the chromosome
		print "Loading Chr$chr start positions into StartHash.\n";
		my %StartHash;  #StartHash{position}[genename1,genename2,ect]
		foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
		{ 
			push($StartHash{$GTFHash{$chrom_name}{$gene_name}[0]},$gene_name);
		}
		
		print "Analyzing Chromosome $chr \n";
		while(<IN>)
		{
			chomp;
			my @InArray = split("\t",$_);
		
			my $startread = $InArray[1];  #start of read
			my $endread = $InArray[2];    #end of read

			foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
			{ 
				if($endread >= $GTFHash{$chrom_name}{$gene_name}[0] && 						$startread <= $GTFHash{$chrom_name}{$gene_name}[1])
				{
					++$GTFHash{$chrom_name}{$gene_name}[3];
				}
			}
		}
		#print GTFHash to Outfile
		foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
		{ 
			my $Chrposition = $chrom_name . ":" . $GTFHash{$chrom_name}{$gene_name}[0] . 					"-" . $GTFHash{$chrom_name}{$gene_name}[1];
			my $genelength = $GTFHash{$chrom_name}{$gene_name}[1] - $GTFHash{$chrom_name}{$gene_name}[0];
			my $RPKMscore = $GTFHash{$chrom_name}{$gene_name}[3] / $genelength;
			$RPKMscore = $RPKMscore / $TotalReads;
			$RPKMscore = sprintf("%.4f", $RPKMscore);		
			print OUT "$gene_name\t$Chrposition\t$genelength\t$GTFHash{$chrom_name}{$gene_name}[2]\t$RPKMscore\n";
		}
	}
	close IN;
	close GTF;
	close OUT;
}

###########################################################################
#                    (15) FPKM from GTF and FPKMWIG                       #
# Scores FPKM (Fragments Per Kilobase per Million mapped reads)           #
# for the given GTF file and FPKMWIG directory/prefix                     #
#                                                                         #
#  Input: 1) Input FPKMWIG prefix                                         #
#         2) Input GTF File                                               #
#         3) Output FPKM Summary Table File                               #
#         4) Column Name for Table                                        #
#                                                                         #
# Output: Returns: FPKM Output File                                       #
###########################################################################

sub fpkm_from_gtf_fpkmwig
{
	# Input
	my ($inputFPKMprefix, $GTFfilename, $OutputName, $ColumnName) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(GTF, "<$GTFfilename") or die "cannot open $GTFfilename GTF infile";
	open(OUT, ">$OutputName") or die "cannot open $OutputName OUT outfile";
	print OUT "Feature_Name\tFull_position\tSub_Positions\tLength\tStrand\t$ColumnName\n";

	my @Chr;             # array that contains all the the names of the mouse chromosomes
	for (my $n = 1; $n< 20; $n++)
	{
		push(@Chr, $n);
	}
	push(@Chr, "M");
	push(@Chr, "X");
	push(@Chr, "Y");

	########################
	# GTF structure create #
	########################

	# Takes a GTF file with lines like this (tab-delimited):
	#     chrom   source    feature   start    end    score   strand
	# ex: chr1  CpG_Island  CpG_100  3204563 3207049    .       +
	# and returns the following data structure:
	#
	# GTF{chrom}->{feature}->{number}->[start,stop,length,strand,summed FPKM over length]
	# 
	# Note: FPKM will be used later in the script but for now will be 0.

	my $source;

	print "Loading GTFHash from GTF file.\n";
	my %GTFHash;
	while(<GTF>)
	{
		chomp;
		my @line = split("\t",$_);
		my $chrom = $line[0];
		$source = $line[1];
		my $feature = $line[2];
		my $start = $line[3];
		my $end = $line[4];
		my $score = $line[5];
		my $strand = $line[6];
	
		# Corrects if start is greater than end
		if($start > $end)
		{
			my $temp = $end;
			$end = $start;
			$start = $temp;
		}
		my $length = $end - $start;
		my $flag = 1;
		while($flag > 0)
		{
			if(exists $GTFHash{$chrom}{$feature}{$flag})
			{
				$flag++;
			}
			else
			{
				push(@{$GTFHash{$chrom}{$feature}{$flag}},$start);
				push(@{$GTFHash{$chrom}{$feature}{$flag}},$end);
				push(@{$GTFHash{$chrom}{$feature}{$flag}},$length);
				push(@{$GTFHash{$chrom}{$feature}{$flag}},$strand);
				push(@{$GTFHash{$chrom}{$feature}{$flag}},0);
				$flag = 0;
			}
		}
	}
	close GTF;

	##############################
	# Calculates and prints FPKM #
	##############################

	while(@Chr)
	{
		# my %FPKMHash;  #holds FPKM from IN file for small selected area
		my $chr = shift(@Chr);

		# opens infile and removes first two lines from it (header lines)
		my $inputfile = $inputFPKMprefix . $chr . ".wig";
		open(IN, "<$inputfile") or die "cannot open $inputfile IN infile";
		<IN>; <IN>;

		my %FeatureHash;
		my $chrom_name = "chr" . $chr;

		print "Loading Chr$chr start positions into FeatureHash.\n";
		my %StartHash;
		# Makes a StartHash and a FeatureHash for the chromosome
		foreach my $feature (keys %{$GTFHash{$chrom_name}}) 
		{ 
			my $feature_start = 999999999999;
			my $feature_end = 0;
			foreach my $segnumber (keys %{$GTFHash{$chrom_name}{$feature}})
			{ 
				if($GTFHash{$chrom_name}{$feature}{$segnumber}[0] < $feature_start)
				{
					$feature_start = $GTFHash{$chrom_name}{$feature}{$segnumber}[0];
				}
				if($GTFHash{$chrom_name}{$feature}{$segnumber}[1] > $feature_end)
				{
					$feature_end = $GTFHash{$chrom_name}{$feature}{$segnumber}[1];
				}		
			}
			# Adds start and end to Feature Hash
			push(@{$FeatureHash{$feature}},$feature_start);
			push(@{$FeatureHash{$feature}},$feature_end);
		
			# Adds start position and corresponding 
			# Note, a single starting position can have multiple features associated with it
			push(@{$StartHash{$feature_start}},$feature);
		}

		my @SortedFeatureArray;
		# Fills @SortedFeatureArray with list (in order of starting position)
		foreach my $startpos (sort { $a <=> $b } keys %StartHash)
		{
			while(exists $StartHash{$startpos}[0])
			{
			push(@SortedFeatureArray,$StartHash{$startpos}[0]);
			shift(@{$StartHash{$startpos}});
			}
		}
	
		undef %StartHash;

		my $lastendpos = 0;

		print "Analyzing Chromosome $chr \n";
		while(<IN>)
		{
			chomp;
			my @InArray = split("\t",$_);
			my $currentpos = $InArray[0];
			my $currentheight = $InArray[1];

			# Ends the loop if there are no more features in this array
			if(exists $SortedFeatureArray[0]){}
			else {last;}
		
			my $sortedfeaturepointer = 0;
			# while current position is not before the start of the current hash
			while($FeatureHash{$SortedFeatureArray[$sortedfeaturepointer]}[0] 					< $currentpos)
			{
				my $feature = $SortedFeatureArray[$sortedfeaturepointer];
			
				# If current position is past the end of the feature,
				# removes feature from SortedFeatureArray
				if($currentpos > $FeatureHash{$feature}[1])
				{
					shift(@SortedFeatureArray);
					# my $feature_finished = shift(@SortedFeatureArray);
					# print "Finished Feature: " , $feature_finished , "\t", $FeatureHash{$feature}[1] , "\n";
					if(exists $SortedFeatureArray[0]){next;}
					else {last;}
				}		
			
				# If current position is before the feature start, exit loop
				if($currentpos < $FeatureHash{$feature}[0])
				{
					last;
				}
		
				# Go through feature area, adding FPKM to areas that overlap with it
				my $segnumber = 1;
				while(exists $GTFHash{$chrom_name}{$feature}{$segnumber})
				{
					if($currentpos > $GTFHash{$chrom_name}{$feature}{$segnumber}[0] && $currentpos < $GTFHash{$chrom_name}{$feature}{$segnumber}[1])
					{
						$GTFHash{$chrom_name}{$feature}{$segnumber}[4] = $GTFHash{$chrom_name}{$feature}{$segnumber}[4] + $currentheight;
					}
					++$segnumber;
				}			
				++$sortedfeaturepointer;
			}
		}
		close IN;
	
		print "Printing $chrom_name to Outfile\n\n";
		#Print FPKM's to Outfile
		foreach my $feature (sort keys %{$GTFHash{$chrom_name}})
		{ 
			my $featurelength = 0;
			my $FPKM = 0;
			my $Chr_total_position = $chrom_name . ":" . $FeatureHash{$feature}[0] . "-" . $FeatureHash{$feature}[1];
			my $Chrpositions = "";
			my $loopcount = 0;
			foreach my $segnumber (keys %{$GTFHash{$chrom_name}{$feature}})
			{
				if ($loopcount > 0)
				{
					$Chrpositions = $Chrpositions . " / ";
				}
				++$loopcount;
				$featurelength = $featurelength + $GTFHash{$chrom_name}{$feature}{$segnumber}[2];
				$FPKM = $FPKM + $GTFHash{$chrom_name}{$feature}{$segnumber}[4];
				$Chrpositions = $chrom_name . ":" . $GTFHash{$chrom_name}{$feature}{$segnumber}[0] . "-" . $GTFHash{$chrom_name}{$feature}{$segnumber}[1];
			}
			if($featurelength > 0)
			{
				$FPKM = (1000 * $FPKM) / $featurelength;
				$FPKM = sprintf("%.4f", $FPKM);	
			} 
			else
			{
				$FPKM = 0;
			}
			# Columns (tab-delimited)
			# Feature  TotalPosition  Positions  FeatureLength(s)  Strand  FPKM
			print OUT "$feature\t$Chr_total_position\t$Chrpositions\t$featurelength\t$GTFHash{$chrom_name}{$feature}{1}[3]\t$FPKM\n";
		}	
	}
	close OUT;

}

1;

