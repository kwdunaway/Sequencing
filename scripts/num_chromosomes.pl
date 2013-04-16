#!/usr/bin/perl
BEGIN {push @INC, "/home/kwdunaway/perl_script";}
use strict; use warnings;
use SeqProcess;

my $infileroot = shift(@ARGV);

 	my @Chr;
	my $filedir = $infileroot;
	$filedir =~ s/^(.*\/)[^\/]*$/$1/; # Gets top directory path
	my @files = glob( $filedir . '*' ); # Gets list of all files in directory
	@files = grep /hr(.+)\.bed/, @files;

	foreach my $file (@files) {
		$file =~ /hr(.+)\.bed/;
		push (@Chr, $1);
	}
	@Chr = sort @Chr;
