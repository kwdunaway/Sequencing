#!/usr/bin/perl 
use strict; use warnings;

die "This script needs the following arguments:
    1) Input
" unless @ARGV == 2;

my $inputfirstbed = shift(@ARGV);
my $inputsecondbed = shift(@ARGV);

my $filedir = $inputfirstbed . "*";
my $firsttotal = `wc -l $filedir | grep "total"`;
$firsttotal =~ s/\D//g;
print "$filedir\n$firsttotal\n";
$filedir = $inputsecondbed . "*";
my $secondtotal = `wc -l $filedir | grep "total"`;
$secondtotal =~ s/\D//g;
print "$filedir\n$secondtotal\n";
my $secondmultiplier = -$firsttotal/$secondtotal;
print "$filedir\n$secondmultiplier\n";
