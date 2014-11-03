#!/usr/bin/perl
# kmer_counter.pl

################################################################################
# 
# Given a raw, de-barcoded Bind-n-Seq reads file,
# counts k-mers using a sliding window of length k
# 
################################################################################

use strict;
use warnings;
use Getopt::Std;

my %opt = (
	'k' => 10,
	'm' => 10,
	'a' => 0,
);

getopts('k:m:a', \%opt);

die "
usage: kmer_counter.pl [options] <raw sequence>
options:
  -k <int> k-mer length [$opt{k}]
  -m <int> minimum number of counts to report [$opt{m}]
  -a       keep ambiguity characters
" unless @ARGV == 1;

my $KLEN = $opt{'k'};
my $MIN = $opt{'m'};

die "unsupported option at this time" if $opt{'a'};

my %kcount;
while (<>) {
	chomp;
	for (my $i = 0; $i < length($_) - $KLEN +1; $i++) {
		my $kmer = substr($_, $i, $KLEN);
		next if $kmer =~ /[^ACGT]/;
		$kcount{$kmer}++;
	}
}

print "kmer_counter output\n";
foreach my $kmer (sort keys %kcount) {
	print  $kmer, "\t", $kcount{$kmer},"\n" if $kcount{$kmer} >= $MIN;
}



__END__



