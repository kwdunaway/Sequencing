#!/usr/bin/perl
use strict; use warnings;

###########################################################################
#               (7) Eliminate Duplicate Reads in BED Files                #
# Checks BED files for duplicate reads and deletes the duplicates         #
#                                                                         #
#  Input: 1) Sorted BED File (with 6 fields)                              #
#         2) Maximum Duplicate Reads allowed (1 for no duplicates)        #
# Output: BED File cleaned of duplicates (replaces input file)            #
#         Returns: 1) Total Number of Duplicate Reads Deleted             #
#                  2) Total Number of Duplicate Read Positions            #
###########################################################################


	# Input
	my ($bedfile, $MaxDupReads) = @ARGV;
	
	##################################################
	#     Global Variables and I/O Initiation        #
	##################################################
	
	my $temp = $bedfile . "_nodups"; # Temporary output file

	open(IN, "<$bedfile") or die "cannot open $bedfile bed file";
	open(OUT, ">$temp") or die "cannot open $temp temp file";

	my $linenum = 1; # Record line number(according to output file)
	# Check if lines are the same; if so, duplicates
	my %data; # Will hold data of every line for checking
	my $DupCount = 1; # Checks for current number of the same read printed to output
	my $TotalDupReads = 0; # Total duplicate reads deleted (exceeded max allowed dups per read)
	my $TotalDupPositions = 0; # Total number of unique reads with duplicates

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
			print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 					$line[4], "\t", $line[5];
			$DupCount = 1;
			$linenum++; 
		}
		elsif ($DupCount < $MaxDupReads)
		{
			print OUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", 					$line[4], "\t", $line[5];
			$DupCount++;
			$linenum++; 
		}
		elsif ($DupCount >= $MaxDupReads)
		{
			
		}
	}

	close(IN);
	close(OUT);

	`rm $bedfile`;
	`mv $temp $bedfile`;

	
		
	

