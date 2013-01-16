#!/usr/bin/perl
use strict; use warnings;

	###########################################################################
	#             (16) Determine Read Length in fastq.gz files                #
	#  Input: Raw file folder (only zipped files and the extension is .fq.gz) #
	# Output: Returns: Filtered and Combined into one .fq file                #
	###########################################################################

	# Input
	my ($rawfqfolder) = @_;

