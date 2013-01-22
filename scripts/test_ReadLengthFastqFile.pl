#!/usr/bin/perl
use strict; use warnings;

###########################################################################
#           (0) Determine Read Length from fastq.gz File                  #
###########################################################################

	# Input
	my ($fastqgzfile) = @ARGV;

	my $readlength = `gunzip -c $fastqgzfile | head -n 2 | tail -n 1 | tr -d '\n'| wc -m | tr -d '\n'`;

	print $readlength, "\n";

# Read count in fastq files

	my $readcount = `gunzip -c $fastqgzfile | wc -l | tr -d '\n'`;
	die "Fastq file format error" if $readcount % 4 != 0;
	print $readcount/4, "\n";
