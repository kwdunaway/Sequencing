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
################################################################################################

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

###########################################################################
#          (0.1) Determine Read Length from Fastq (.fq) File              #
###########################################################################

sub fastq_readlength 
{
	my ($fastqfile) = @_; # Input is "filename.fq"
	print "Read Length from .fq File: Obtaining read length from " , $fastqfile , "\n";
	# Determines read length from the second line of the file
	my $readlength = `head -n 2 $fastqfile | tail -n 1 | tr -d '\n'| wc -m | tr -d '\n'`;	
	print "Read Length from .fq File: Finished. Read length is ", $readlength ," \n\n";
	return $readlength;
}

###########################################################################
#            (0.2) Determine Read Count from Fastq (.fq) File             #
###########################################################################

sub fastq_readcount
{
	my ($fastqfile) = @_; #Input is "filename.fq"
	print "Read Count from .fq File: Obtaining read count from " , $fastqfile , "\n";
	# Uses word count to find the read count and divides by 4 to return the actual read count
	my $readcount = `wc -l $fastqfile | tr -d '\n'`;
	die "Error: Read Count from .fq File: Fastq file format error" if $readcount % 4 != 0;
	my $truereadcount = $readcount/4;
	print "Read Count from .fq File: Finished. Read count is " , $truereadcount , "\n";
	return $truereadcount;
}

###########################################################################
#           (0.3) Determine Read Count from Eland Extended File           #
###########################################################################

sub elandext_readcount
{
	my ($elandextfile) = @_; #Input is Eland Extended file name
	print "Read Count from Eland Extended File: Obtaining read count from " , $elandextfile , "\n";
	# Uses word count to return the read count
	my $readcount = `wc -l $elandextfile`;
	print "Read Count from Eland Extended File: Finished. Read count is " , $readcount , "\n";
	return $readcount;
}

###########################################################################
#                     (1) Header and Path Modifier                        #
###########################################################################

sub add_path 
{
	my ($addtoPATH) = @_; # Input is path to add
	print "Add to PATH: Adding " , $addtoPATH , " to PATH for rest of this script.\n";
	$ENV{'PATH'} = $ENV{'PATH'} . ":" . $addtoPATH;
	print "Add to PATH: Finished. Added " , $addtoPATH , " to PATH for rest of this script.\n\n";	
}

###########################################################################
#              (2) Combine and Filter Zipped Fastq Files                  #
#  Input: Raw data file folder (extension is .fq.gz and only .fq.gz)      #
# Output: Returns: Filtered and Combined into one .fq file                #
###########################################################################

sub filter_zip 
{
	my ($rawfqfolder) = @_; # Input is path to folder containing data with extension .fq.gz
	my $filtered_fastq = $rawfqfolder . "filtered.fq"; # Filtered file is named "filtered.fq" in input folder
	print "Filter zipped fastq files: Filtering $rawfqfolder files and outputting to $filtered_fastq\n";
	# Uses gunzip and grep to filter the data into "filtered.fq"
	my $comline = "gunzip -c " . $rawfqfolder . "*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v \"^--\$\" >  " . $filtered_fastq;
	`$comline`;
	print "Filter zipped fastq files: Finished. Output to $filtered_fastq\n\n";
#	return $filtered_fastq;
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
	# Input (see subroutine description)
	my ($ExperimentTopDir, $BowtiePrefix, $mm9path, $filtered_fastq) = @_;

	# Making folder to contain output
	print "Bowtie: Making $ExperimentTopDir directory in current directory\n";
	`mkdir $ExperimentTopDir\n`;

	# File names for non-aligned and aligned reads files (e.g. $ExperimentTopDir/$BowtiePrefix_NonAligned.fq)
	my $nonalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_NonAligned.fq";
	my $alignedpreseparationfile = $ExperimentTopDir . $BowtiePrefix . "_alignedpreseparation.txt";

	# Run Bowtie
	print "Bowtie: Running Bowtie, separating aligned reads to $alignedpreseparationfile and non-aligned reads to $nonalignedreadsfile\n";
	`bowtie -p 4 -M 1 -k 1 --chunkmbs 256 --strata --best --un $nonalignedreadsfile $mm9path $filtered_fastq $alignedpreseparationfile\n\n`;
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
	# Input (see subroutine description)
	my ($ExperimentTopDir, $BowtiePrefix, $alignedpreseparationfile) = @_;

	# File names for unique and repeat reads files (e.g. $ExperimentTopDir/$BowtiePrefix_Uniq.txt)
	my $uniqalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Uniq.txt";
	my $repalignedreadsfile = $ExperimentTopDir . $BowtiePrefix . "_Repeat.txt";

	print "Separate Repeats and Uniques: separating aligned reads between unique reads to $uniqalignedreadsfile and repeat reads to $repalignedreadsfile\n";

	# Open Files
	open(IN, "<$alignedpreseparationfile") or die "Error: Separate Repeats and Uniques: cannot open $alignedpreseparationfile infile";
	open(UNIQOUT, ">$uniqalignedreadsfile") or die "Error: Separate Repeats and Uniques: cannot open $uniqalignedreadsfile outfile";
	open(REPOUT, ">$repalignedreadsfile") or die "Error: Separate Repeats and Uniques: cannot open $repalignedreadsfile outfile";

	while (<IN>) {
		chomp;
		my @line = split ("\t", $_);
		# Check if column 7 is greater than 0, yes means repeat read, no means unique read
		if($line[6] > 0) {
			print REPOUT $_, "\n";
		}
		else {
			print UNIQOUT $_, "\n";
		}
	}

	# Close Files
	close IN;
	close UNIQOUT;
	close REPOUT;

	print "Separate Repeats and Uniques: Finished.\n";
}

###########################################################################
#                   (5) Eland Extended Format to BED                      #
#   *             Note: subroutine calls the "Sort BED"               *   #
#   *        and "Eliminate Duplicate BED Reads" subroutines          *   #
#                                                                         #
#  This subroutine takes sequences in Eland Extended format (also known   #
#  as s_#_export.txt files) and produces multiple BED files which can     #
#  be used to analyze the data. Requires position in array of chromosome  #
#  number, position, and strand (+,-)                                     #
#                                                                         #
#  Input: 1) Input file name                                              #
#         2) Top folder path of experiment (Folder to contain)            #
#         3) Output directory name                                        #
#         4) Starting Read Length (Current)                               #
#         5) Final Read Length (Final to replace current)                 #
#         6) Chromosome array number                                      #
#         7) Position array number                                        #
#         8) Strand array number                                          #
#         9) Maximum Duplicate Reads (1 for no duplicates)                #
#                                                                         #
# Output: Creates directory containing:                                   #
#           1) BED files for each chromosome                              #
#           2) Statistics on all reads                                    #
#           3) Gzipped non-unique and non-mappable reads files            #
#           4) Unknown reads file                                         #
###########################################################################

sub elandext_to_bed 
{
	print "\nElandExt to BED: Beginning conversion of Eland Extended format to BED format\n";
	# Input (see subroutine description)
	my ($infile, $ExperimentTopDir, $FilePrefix, $basereadlength, $finalreadlength, $chr, $pos, $strand, $MaxDupReads) = @_;

	# Variable Declarations
	my $outdir = $ExperimentTopDir . $FilePrefix . "_bed"; # Output directory name
	my $minusstrandlength = $finalreadlength - $basereadlength; # Calculation of the minus strand length
	my %Count; # Keys = chromosome number, Values = number of mapped reads for specific chromosome
	my %Files; # Contains files with each chromosome
	my $totalcount = 0; # Total number of mapped reads

	# Makes Output Directory if does not exist
	print "ElandExt to BED: Making $outdir directory if does not exist\n";
	if (! -d $outdir) 
	{ 
		`mkdir $outdir`;
		if (! -d $outdir) { die "Error: ElandExt to BED: directory ($outdir) still does not exist"} 
	}

	open(IN, "<$infile") or die "Error: ElandExt to BED: cannot open $infile infile"; 
	print "ElandExt to BED: Processing data into each chromosome output file in BED format\n";

	############################################################################
	#                      Processing Eland Extended File                      #
	############################################################################
	while (<IN>)
	{
		chomp;
		my @array = split("\t", $_); # Splitting tab-delimited data into array
		my $chrom = $array[$chr];
		my $readstrand = $array[$strand];
		$totalcount++; # Add to total mapped reads (each line is a read)

		#If output bed file has not been created for the chromosome, create it
		if(! exists $Count{$chrom}){
			$Count{$chrom} = 0; # Initialize read count for this chromosome 
			my $filename = $outdir . "/" . $FilePrefix . "_" . $chrom . ".bed";
			open($Files{$chrom}, ">$filename") or die "cannot open $filename chromosome output file";	
		}
		my $printfile = $Files{$chrom};
		$Count{$chrom}++; # Add to number of mapped reads for the chromosome

		# Print line to bed file according to whether the strand is + or -
		if($readstrand eq "+"){
			print {$Files{$chrom}} $chrom , "\t" , $array[$pos] , "\t" , $array[$pos] + $finalreadlength, "\t", $FilePrefix , "\t", "0", "\t" , $readstrand , "\n";
		}
		elsif($readstrand eq "-"){
			print {$Files{$chrom}} $chrom , "\t" , $array[$pos] - $minusstrandlength, "\t" , $array[$pos] + $basereadlength, "\t", $FilePrefix , "\t", "0", "\t" , $readstrand , "\n";
		}
		else {die "Error: ElandExt to BED: Strand in read is not + nor -, it is $readstrand";}
	}
	close IN;
	print "ElandExt to BED: Finished outputting data to each BED file\n";


	############################################################################
	#                  Printing statistics to Stats Outfile                    #
	############################################################################
	
	# Statistics file name is Stats_fileprefix.txt in the specified folder
	my $filename = $ExperimentTopDir . "/" . "Stats_" . $FilePrefix . ".txt";
	print "ElandExt to BED: Printing statistics to ", $filename, "\n";
	open(STATS, ">$filename") or die "Error: ElandExt to BED: cannot open $filename statistics output file";	
	# While printing reads mapped to each chromosome, sort each bed file and eliminate duplicate reads
	foreach my $key (sort keys %Files) {
		print STATS "Number of reads mapped " , $key , " is:\t" , $Count{$key} , "\n";
		close $Files{$key};
		my $bedfile = $outdir . "/" . $FilePrefix . "_" . $key . ".bed"; # Name of BED file
		sort_bed($bedfile); # Sort BED file
		eliminate_bed_dups($bedfile, $MaxDupReads); # Eliminate unwanted duplicate reads
	}
	print STATS "\nTotal Number of Total mapped reads is:\t", $totalcount, "\n";
	close STATS;
	print "ElandExt to BED: Finished conversion of Eland Extended to BED format\n\n";
}

###########################################################################
#                           (6) Sort BED File                             #
#  Input: Unsorted BED file                                               #
# Output: Sorted BED file (replaces unsorted file)                        #
###########################################################################

sub sort_bed
{
	my ($bedfile) = @_;
	my $temp = $bedfile . "_sorted"; # Creates temporary name for the sorted bed file
	print "Sort BED: Sorting BED file $bedfile\n";
	`sort -n +1 -2 $bedfile > $temp`; # Numerically sort bed file contents using command line
	# Replace unsorted bed file with sorted bed file
	`rm $bedfile`; 
	`mv $temp $bedfile`; 
	print "Sort BED: Finished.\n";
}

###########################################################################
#               (7) Eliminate Duplicate Reads in BED Files                #
# Checks BED files for duplicate reads and deletes the duplicates.        #
# BED file must be sorted.                                                #
#                                                                         #
#  Input: 1) Sorted BED File (with 6 fields)                              #
#         2) Maximum Duplicate Reads allowed (1 for no duplicates)        #
# Output: BED File cleaned of duplicates (replaces input file)            #
#         Returns: 1) Total Number of Duplicate Reads Deleted             #
#                  2) Total Number of Duplicate Read Positions            #
###########################################################################

sub eliminate_bed_dups
{
	# Input (see subroutine description)
	my ($bedfile, $MaxDupReads) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################
	
	my $temp = $bedfile . "_nodups"; # Temporary name for output file

	open(IN, "<$bedfile") or die "Error: Eliminate Duplicate Reads: cannot open $bedfile bed file";
	open(OUT, ">$temp") or die "Error: Eliminate Duplicate Reads: cannot open $temp temp file";

	my $linenum = 1; # Record line number(according to output file)
	# Check if lines are the same; if so, duplicates
	my %data; # Will hold data of every line for checking
	my $DupCount = 1; #Checks for current number of the same read printed to output
	my $TotalDupReads = 0; # Total duplicate reads deleted (exceeded max allowed dups per read)
	my $TotalDupPositions = 0; # Total number of unique reads with duplicates

	##################################################
	#            Checking for Duplicates             #
	##################################################

	print "Eliminate Duplicate Reads: Eliminating duplicate reads in $bedfile\n";
	# Retrieve data for first line and split into array
	$_ = <IN>;
	$data{$linenum} = $_;
	my @line = split("\t", $_);
	print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5];
	$linenum++;

	while (<IN>)
	{
		$data{$linenum} = $_;
		my @line = split("\t", $_);

		# If current line does not equal previous line, print to output
		if ($data{$linenum-1} ne $data{$linenum})
		{
			print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 				$line[4], "\t", $line[5];
			$DupCount = 1; # Resets duplicate counter
			$linenum++; 
		}
		# If current line equals previous line but maximum allowed duplicates is not reached, print to output
		elsif ($DupCount < $MaxDupReads)
		{
			print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 					$line[4], "\t", $line[5];
			$DupCount++;
			$linenum++; 
		}
		# If current line equals previous line and maximum allowed duplicates has been reached, do not print to output and do not advance line number. Add to data on number of duplicate reads deleted.
		elsif ($DupCount >= $MaxDupReads)
		{
			$TotalDupReads++; # Total duplicate reads deleted
			if ($DupCount == $MaxDupReads) # Count once per position
			{
				$TotalDupPositions++; # Total number of unique reads with duplicates
			}
			$DupCount++;
		}
	}

	close(IN);
	close(OUT);

	# Replace input file with new file
	`rm $bedfile`;
	`mv $temp $bedfile`; 
	
	print "Eliminate Duplicate Reads: Finished. Total duplicate reads deleted is $TotalDupReads. Total number of unique reads with duplicates is $TotalDupPositions.\n";

	# Returns variables with statistics on reads deleted
	return ($TotalDupReads, $TotalDupPositions);
}

###########################################################################
#                   (8) Extend Read Length of BED Files                   #
# Takes a folder of bed files and creates a new one with read length      #
# extended directionally based on arguments input on command line.        #
# Input files must end with "chr*.bed" (* is chromosome number)           #
#                                                                         #
#  Input: 1) Input Bed File prefix (ex: DY specifies DY_chr*.bed)         #
#         2) Output Bed prefix (ex: NewDY will make NewDY_chr*.bed)       #
#         3) Read Length extension (ex: 97)(added to old read length)     #
#                                                                         #
# Output: Extended Read Length BED Files                                  #
###########################################################################

sub extend_bed_read_length
{
	# Input (see subroutine description)
	my ($inputprefix, $outputbedprefix, $readlengthextension) = @_;

	print "Extend Read Length: Beginning extension of read lengths in BED files by $readlengthextension\n";

	die "Error: Extend Read Length: Input BED prefix ($inputprefix) is same as the output BED prefix ($outputbedprefix).\n" if $inputprefix eq $outputbedprefix;

	# Scan BED directory for number of chromosome files
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	# For every chromosome number detected, look for an input bed file with specified chromosome number and extend read length
	while(@Chr)
	{
		# File name creation and opening of files
		my $chr = shift(@Chr);
		my $inputfile = $inputprefix . "_chr" . $chr . ".bed";
		open(IN, "<$inputfile") or die "Error: Extend Read Length: cannot open $inputfile infile";
		my $outfile = $outputbedprefix . "_chr" . $chr . ".bed";
		open(OUT, ">$outfile") or die "Error: Extend Read Length: cannot open $outfile outfile";

		print "Extend Read Length: Extending read lengths in $inputfile\n";

		# Open input bed file pertaining to the specified chromosome number
		while(<IN>)
		{
			chomp;
    			my @line = split ("\t", $_);
			my $start = $line[1]; # Column 2: Start location of input read
			my $end = $line[2]; # Column 3: End location of input read
			# Extend according to if strand is + or -
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
			else {die "Error: Extend Read Length: $line[5] does not equal + or - \n";}
		}

		close(IN);
		close(OUT);
	}
	print "Extend Read Length: Finished read length extensions\n";
}

###########################################################################
#            (9) Change Read Length of BED Files (Choose Read Length)     #
# Takes a folder of bed files and creates a new one with read length      #
# based on arguments input on command line. This version of the "Extend   #
# Read Length" subroutine allows choosing the final read length as        #
# opposed to extending from the base.                                     #
# Input files must end with "chr*.bed" (* is chromosome number)           #
#                                                                         #
#  Input: 1) Input Bed File prefix (ex: DY specifies DY_chr*.bed)         #
#         2) Output Bed prefix (ex: NewDY will make NewDY_chr*.bed)       #
#         3) New Read Length (ex: 97)                                     #
#                                                                         #
# Output: New Read Length BED Files                                       #
###########################################################################

sub change_bed_read_length
{
	# Input (see subroutine description)
	my ($inputprefix, $outputbedprefix, $readlengthextension) = @_;

	print "Change Read Length: Changing read lengths in BED files to $readlengthextension\n";

	die "Error: Change Read Length: Input BED prefix ($inputprefix) is same as the output BED prefix ($outputbedprefix).\n" if $inputprefix eq $outputbedprefix;

	# Scan BED directory for number of chromosome files
	my @Chr;
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	# For every chromosome number detected, look for an input bed file with specified chromosome number and change read length
	while(@Chr)
	{
		# File name creation and opening of files
		my $chr = shift(@Chr);
		my $inputfile = $inputprefix . "_chr" . $chr . ".bed";
		open(IN, "<$inputfile") or do {next;};
		my $outfile = $outputbedprefix . "_chr" . $chr . ".bed";
		open(OUT, ">$outfile") or die "Error: Change Read Length: cannot open $outfile outfile";

		print "Change Read Length: Changing read lengths in $inputfile\n";

		# Open input bed file pertaining to the specified chromosome number
		while(<IN>)
		{
			chomp;
    			my @line = split ("\t", $_);
			my $start = $line[1]; # Column 2: Start location of input read
			my $end = $line[2]; # Column 3: End location of input read
			# Change read length according to if strand is + or -
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
			else {die "Error: Change Read Length: $line[5] does not equal + or - \n";}
		}

		close(IN);
		close(OUT);
	}
	print "Change Read Length: Finished read length change\n";
}

###########################################################################
#              (10) BED Directory to Variable Step WIG                    #
# Converts a directory of BED files to a directory of WIG files           #
# Input files must end with "chr*.bed" (* is chromosome number)           #
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
	# Input (see subroutine description)
	my ($inputprefix, $outfileroot, $wignameroot, $color) = @_;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	# Scan BED directory for number of chromosome files
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	####################################################################################
	# Grabs the information from each line and assigns it to the appropriate variables #
	####################################################################################

	print "\n\nBED Directory to Variable Step WIG: Starting Bed to Wig conversion of files with prefix $inputprefix:\n";

	# Run for every chromosome number detected
	while(@Chr)
	{
		# File name creation and opening of files
		my $infile = $inputprefix . "_chr" . $Chr[0] . ".bed";
		open(IN, "<$infile") or die "Error: BED Directory to Variable Step WIG: cannot open $infile infile"; #opens input file to be read (must be .bed)
		my @inputprefix = split(".bed", $infile);
		my $outfile = $outfileroot . "_chr" . $Chr[0] . ".wig";
		open(OUT, ">$outfile") or die "Error: BED Directory to Variable Step WIG: cannot open $outfile outfile"; #opens output file to write to (.wig)

		# Prints the head of the track (necessary for genome browser to read file properly) 			(customizable through terminal)
		my $wigname = $wignameroot . "_chr" . $Chr[0];
		print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", 
			$wigname, "\" description=\"", $wigname, "\" color=", $color, "\n";

		my %PosVal;              # hash that contains the peak heights, kept small
		my $position = 0;        # position of window
		my $chromcheck;          # variable that checks each chromosome to make sure they 							are the same chromosome
		my @line;                # temp array used to retrieve the information of each line
		my $chrom;               # chromosome of current line
		my $startread = 0;       # start of read at current line
		my $endread = 0;         # end of read at current line 						(does not include this position)

		# Open input bed file pertaining to the specified chromosome number
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
				die "Error: BED Directory to Variable Step WIG: You have different chromosomes in this bed file (Program ended before completion)/n" unless $chrom eq $chromcheck;
	    
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
		print "BED Directory to Variable Step WIG: Finished with Chromosome", $Chr[0] ,"\n";
		shift(@Chr);
	}

}

###########################################################################
#                  (11) Variable Step WIG to FPKM WIG                     #
# Calculates FPKM (Fragments Per Kilobase per Million mapped reads)       #
# from raw WIG files using read length and read count (in millions)       #
# Input files must end with "chr*.wig" (* is chromosome number)           #
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
	# Input (see subroutine description)
	my ($inputprefix, $outprefix, $readlength, $readcount) = @_;

	print "Variable Step WIG to FPKM WIG: Converting $inputprefix Variable Step WIG files to $outprefix FPKM WIG files\n";

	die "Error: Variable Step WIG to FPKM WIG: Input WIG prefix ($inputprefix) is same as the output WIG prefix ($outprefix).\n" if $inputprefix eq $outprefix;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	# Scan WIG directory for number of chromosomes
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.wig/, @files; # Takes only the files with "hr*.wig"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.wig/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	#############################################
	# Calculate and print OUT FPKM from file(s) #
	#############################################

	# Run for every chromosome number detected
	while(@Chr)
	{
		# File name creation and opening of files
		my $chr = shift(@Chr);
		print "Variable Step WIG to FPKM WIG: Now Converting: Chr$chr\n";
		my $inputWIG = $inputprefix . "_chr" . $chr . ".wig";
		my $outfile = $outprefix . "_chr" . $chr . ".wig";
		open(INWIG, "<$inputWIG") or die "Error: Variable Step WIG to FPKM WIG: cannot open $inputWIG infile";
		open(OUT, ">$outfile") or die "Error: Variable Step WIG to FPKM WIG: cannot open $outfile outfile";

		# Print head of track (first two lines) from input wig to output wig
		my $lin = <INWIG>;
		print OUT $lin;
		$lin = <INWIG>;
		print OUT $lin;

		# For rest of the input file, convert to FPKM and print to output wig
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
	print "Variable Step WIG to FPKM WIG: Completed conversion of $inputprefix Variable Step WIG to $outprefix FPKM WIG\n";
}

###########################################################################
#                    (12) BED Directory to FPKM WIG                       #
#   * Note: subroutine calls the "BED Directory to Variable Step WIG" *   #
#   *        and "Variable Step WIG to FPKM WIG" subroutines          *   #
#                                                                         #
# Converts a directory of BED files to a directory of FPKM WIG files.     #
# Input files must end with "chr*.bed" (* is chromosome number)           #
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
	# Input (see subroutine description)
	my ($inputprefix, $outfileroot, $wignameroot, $color, $readlength, $readcount) = @_;


	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	# Scan BED directory for number of chromosome files
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	####################################################################################
	# Grabs the information from each line and assigns it to the appropriate variables #
	####################################################################################

	print "\n\nBED Directory to FPKM WIG: Starting Bed to Wig conversion of files with prefix $inputprefix:\n";

	# Run for every chromosome number detected
	while(@Chr)
	{
		# File name creation and opening of files
		my $infile = $inputprefix . "_chr" . $Chr[0] . ".bed";
		open(IN, "<$infile") or die "Error: BED Directory to FPKM WIG: cannot open $infile infile"; #opens input file to be read (must be .bed)
		my @inputprefix = split(".bed", $infile);
		my $outfile = $outfileroot . "_chr" . $Chr[0] . ".wig";
		open(OUT, ">$outfile") or die "Error: BED Directory to FPKM WIG: cannot open $outfile outfile"; #opens output file to write to (.wig)

		# Prints the head of the track (necessary for genome browser to read file properly) 			(customizable through terminal)
		my $wigname = $wignameroot . "_chr" . $Chr[0];
		print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", 
			$wigname, "\" description=\"", $wigname, "\" color=", $color, "\n";

		my %PosVal;              # hash that contains the peak heights, kept small
		my $position = 0;        # position of window
		my $chromcheck;          # variable that checks each chromosome to make sure they 							are the same chromosome
		my @line;                # temp array used to retrieve the information of each line
		my $chrom;               # chromosome of current line
		my $startread = 0;       # start of read at current line
		my $endread = 0;         # end of read at current line 						(does not include this position)

		# Open input bed file pertaining to the specified chromosome number
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
				die "Error: BED Directory to FPKM WIG: You have different chromosomes in this bed file (Program ended before completion)/n" unless $chrom eq $chromcheck;
	    
				# adds height to PosVal for the sequence found
				my $addcount = $startread;
				while ($addcount < $endread)
				{
					if (exists $PosVal{$addcount})
					{
						$PosVal{$addcount} = $PosVal{$addcount} + 1000 / ($readlength * $readcount);
						$PosVal{$addcount} = sprintf("%.4f",$PosVal{$addcount});
					}
					else
					{
						$PosVal{$addcount} = 1000 / ($readlength * $readcount);
						$PosVal{$addcount} = sprintf("%.4f",$PosVal{$addcount});
					}
					$addcount = $addcount + 1;
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
		print "BED Directory to FPKM WIG: Finished with Chromosome", $Chr[0] ,"\n";
		shift(@Chr);
	}
}

###########################################################################
#                         (13) Visualize FPKM WIG                         #
# Combines FPKM (Fragments Per Kilobase per Million mapped reads)         #
# files given multipliers for them.                                       #
# Input files must end with "chr*.wig" (* is chromosome number)           #
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
	# Input (see subroutine description)
	my ($inputprefix, $outprefix, $stepsize, $color, $tracknameprefix) = @_;

	print "Visualize FPKM WIG: Visualizing $inputprefix FPKM WIG files\n";

	die "Error: Visualize FPKM WIG: Input WIG prefix ($inputprefix) is same as the output WIG prefix ($outprefix).\n" if $inputprefix eq $outprefix;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	# Scan WIG directory for number of chromosomes
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.wig/, @files; # Takes only the files with "hr*.wig"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.wig/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

	#############################################
	# Visualize and print OUT FPKM from file(s) #
	#############################################

	# Run for every chromosome number detected
	while(@Chr)
	{
		# File name creation and opening of files
		my $chr = shift(@Chr);
		print "Visualize FPKM WIG: Now Visualizing: Chr$chr\n";
		my $inputWIG = $inputprefix . "_chr" . $chr . ".wig";
		open(IN, "<$inputWIG") or die "Error: Visualize FPKM WIG: cannot open $inputWIG IN infile";
		# Delete first two header lines
		my $lin = <IN>;
		$lin = <IN>;
		my $outfile = $outprefix . "_chr" . $chr . ".wig";	
		open(OUT, ">$outfile") or die "Error: Visualize FPKM WIG: cannot open $outfile OUT outfile";

		# Create head of track in output file
		my $trackname = $tracknameprefix . "_" . $chr;
		print OUT "track type=wiggle_0 visibility=full autoScale=off name=\"", $trackname, "\" description=\"", $trackname, "\" color=", $color, "\n";
		print OUT "variableStep chrom=chr", $chr," span=", $stepsize ,"\n";

		my %HeightHash; # Keys = position, Values = height
		$lin = <IN>; # Reads first line of input wig
		my @firstline = split("\t",substr($lin, 0, -1));
		my $currentpos = $firstline[0]; # Initializes current position at first position
		my $startpos = $firstline[0]; # Initializes start position at first position
		$HeightHash{$firstline[0]}=$firstline[1]; # Adds first position and height to hash

		while(<IN>)
		{
			chomp;
			my @line = split("\t",$_);
			$HeightHash{$line[0]}=$line[1]; # Adds position and height to hash
			$currentpos = $line[0]; # Current position at current line
			if ($startpos + $stepsize <= $currentpos) # When gap between start position and current position is large enough (of at least size $stepsize), calculate height and print
			{
				my $height=0; # Initialize height for window
				# For every base in the window, add each height to the total if exists
				for(my $n = $startpos; $n < ($startpos + $stepsize); $n++)
				{
					if(exists $HeightHash{$n}) # May not exist if no reads at that base
					{	
						$height = $height + $HeightHash{$n};
						delete $HeightHash{$n};
					}
				}
				$height = $height / $stepsize; # Average the height over the window size
				$height = sprintf("%.2f", $height);
				print OUT $startpos , "\t" , $height , "\n"; # Print start position and height to output
				$startpos = $currentpos; # New start position from current position
			}
		}

		# Run a last window to check if any bases were missed from the last position of the file
		my $height=0;
		for(my $n = $startpos; $n < ($startpos + $stepsize); $n++)
		{
			if(exists $HeightHash{$n})
			{	
				$height = $height + $HeightHash{$n};
			}
		}
		$height = $height / $stepsize;
		$height = sprintf("%.2f", $height);
		print OUT $startpos , "\t" , $height , "\n";

		close IN;
		close OUT;
	
		print "Visualize FPKM WIG: Finished making $outfile, now zipping it.\n";
		my $commandline = "gzip $outfile";
		`$commandline`;
	}
	print "Visualize FPKM WIG: Finished Visualizing $inputprefix FPKM WIG files\n";
}

###########################################################################
#                       (14) RPKM From BED and GTF                        #
# Scores RPKM (Reads Per Kilobase per Million mapped reads)               #
# for the given GTF file and Bed directory/prefix.                        #
# Input bed files must end with "chr*.bed" (* is chromosome number)       #
#                                                                         #
#  Input: 1) Input BED prefix						  #
#         2) Input GTF file                                               #
#         3) Output RPKM table file name                                  #
#         4) Total Reads (in millions)                                    #
#                                                                         #
# Output: RPKM Output File                                                #
###########################################################################

sub rpkm_from_bed
{
	print "RPKM From BED and GTF: Loading GTFHash from GTF file\n";
	# Input (see subroutine description)
	my ($inputprefix, $GTFfilename, $OutputName, $TotalReads) = @_;

	print "RPKM From BED and GTF: Scoring RPKM from $inputprefix BED files and $GTFfilename GTF file\n";

	die "Error: RPKM From BED and GTF: Input GTF file ($GTFfilename) is same as the output RPKM file ($OutputName).\n" if $GTFfilename eq $OutputName;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(GTF, "<$GTFfilename") or die "Error: RPKM From BED and GTF: cannot open $GTFfilename GTF infile";
	open(OUT, ">$OutputName") or die "Error: RPKM From BED and GTF: cannot open $OutputName outfile";
	print OUT "Gene_Name\tChrom_pos\tGene_length\tStrand\tRPKM\n";

	# Scan BED directory for number of chromosome files
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files; # Takes only the files with "hr*.bed"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

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

	print "RPKM From BED and GTF: Loading GTFHash from GTF file\n";
	my %GTFHash; # GTF{chrom}->{gene_name}->[start,stop,strand,count]

	while(<GTF>)
	{
		chomp;
		# Splitting data into array and variables
		my @line = split("\t",$_);
		my $chrom = $line[0];
		my $source = $line[1];
		my $feature = $line[2];
		my $start = $line[3];
		my $end = $line[4];
		my $score = $line[5];
		my $strand = $line[6];
		my $frame = $line[7];
		# Splitting attribute info into array and variables
		my @attribute = split("\"",$line[8]);
		my $gene_id = $attribute[1];
		my $gene_name = $attribute[3];
		my $p_id = $attribute[5];
		my $transcript_id = $attribute[7];
		my $tss_id = $attribute[9];	
		if($feature eq "exon") # Only add exons to GTFHash
		{
			# If gene already exists in hash, use the smallest position as the start and the largest position as the end
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
			# Otherwise, gene does not exist in hash yet, so add gene information to hash
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

	# Run for every chromosome number detected in bed directory
	while(@Chr)
	{
		my $chr = shift(@Chr);
		my $chrom_name = "chr" . $chr;

		# opens infile
		my $inputfile = $inputprefix . "_chr" . $chr . ".bed";
		open(IN, "<$inputfile") or die "Error: RPKM From BED and GTF: cannot open $inputfile IN infile";

		# Makes a StartHash for the chromosome
		print "RPKM From BED and GTF: Loading Chr$chr start positions into StartHash\n";
		my %StartHash;  # StartHash{position}[genename1,genename2,ect]
		# Add all genes at that start position to the position in the hash
		foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
		{ 
			push(@{$StartHash{$GTFHash{$chrom_name}{$gene_name}[0]}},$gene_name);
		}
		
		print "RPKM From BED and GTF: Analyzing Chromosome $chr \n";
		while(<IN>)
		{
			chomp;
			my @InArray = split("\t",$_);
		
			my $startread = $InArray[1];  #start of read
			my $endread = $InArray[2];    #end of read

			# Run for every gene in the hash
			foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
			{ 
				# If read maps to any part of the gene, add to count
				if($endread >= $GTFHash{$chrom_name}{$gene_name}[0] && 						$startread <= $GTFHash{$chrom_name}{$gene_name}[1])
				{
					++$GTFHash{$chrom_name}{$gene_name}[3];
				}
			}
		}
		# Print GTFHash to Outfile, for every gene
		foreach my $gene_name (keys %{$GTFHash{$chrom_name}})
		{ 
			# $Chrposition = chr*:start-end
			my $Chrposition = $chrom_name . ":" . $GTFHash{$chrom_name}{$gene_name}[0] . 					"-" . $GTFHash{$chrom_name}{$gene_name}[1];
			# Length of the gene (End position - start position)
			my $genelength = $GTFHash{$chrom_name}{$gene_name}[1] - $GTFHash{$chrom_name}{$gene_name}[0];
			# RPKM score = (Count/Gene Length) / Total Reads (in Millions)
			my $RPKMscore = $GTFHash{$chrom_name}{$gene_name}[3] / $genelength;
			$RPKMscore = $RPKMscore / $TotalReads;
			$RPKMscore = sprintf("%.4f", $RPKMscore);
			# Print to output (tab-delimited)
			# 1		2		3		4	5
			# gene name	chr*:start-end	gene length	strand	RPKM score		
			print OUT "$gene_name\t$Chrposition\t$genelength\t$GTFHash{$chrom_name}{$gene_name}[2]\t$RPKMscore\n";
		}
	}
	close IN;
	close GTF;
	close OUT;
	print "RPKM From BED and GTF: Finished with $OutputName RPKM file\n";
}

###########################################################################
#                    (15) FPKM from GTF and FPKMWIG                       #
# Scores FPKM (Fragments Per Kilobase per Million mapped reads)           #
# for the given GTF file and FPKMWIG directory/prefix                     #
# Input files must end with "chr*.wig" (* is chromosome number)           #
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
	# Input (see subroutine description)
	my ($inputprefix, $GTFfilename, $OutputName, $ColumnName, @Chromosomes) = @_;

	print "FPKM from GTF and FPKMWIG: Scoring FPKM from $inputprefix FPKM WIG files and $GTFfilename GTF file\n";

	die "Error: FPKM from GTF and FPKMWIG: Input GTF file ($GTFfilename) is same as the output FPKM file ($OutputName).\n" if $GTFfilename eq $OutputName;

	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################

	open(GTF, "<$GTFfilename") or die "Error: FPKM from GTF and FPKMWIG: cannot open $GTFfilename GTF infile";
	open(OUT, ">$OutputName") or die "Error: FPKM from GTF and FPKMWIG: cannot open $OutputName OUT outfile";
	print OUT "Feature_Name\tFull_position\tSub_Positions\tLength\tStrand\t$ColumnName\n";

	# Scan WIG directory for number of chromosomes
	my @Chr; # Holds all the chromosome numbers (e.g. 19, M)
	my $filedir = $inputprefix;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path from the input prefix
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.wig/, @files; # Takes only the files with "hr*.wig"

	foreach my $file (@files) {
		$file =~ /hr(.+)\.wig/; # For each of those files, extract the chromosome number
		push (@Chr, $1); # Add to list of chromosome numbers
	}
	@Chr = sort @Chr; # Sort list of chromosome numbers

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

	print "FPKM from GTF and FPKMWIG: Loading GTFHash from GTF file.\n";
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
		my $inputfile = $inputprefix . $chr . ".wig";
		open(IN, "<$inputfile") or die "Error: FPKM from GTF and FPKMWIG: cannot open $inputfile IN infile";
		<IN>; <IN>;

		my %FeatureHash;
		my $chrom_name = "chr" . $chr;

		print "FPKM from GTF and FPKMWIG: Loading Chr$chr start positions into FeatureHash.\n";
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

		print "FPKM from GTF and FPKMWIG: Analyzing Chromosome $chr \n";
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
	
		print "FPKM from GTF and FPKMWIG: Printing $chrom_name to Outfile\n\n";
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
	print "FPKM from GTF and FPKMWIG: Finished with $OutputName FPKM file\n";
}

1;

