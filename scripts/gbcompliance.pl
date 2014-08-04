#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Roy Chu and Keith Dunaway
# Email: rgchu@ucdavis.edu
# Date: 8-4-2014
# Script Name: gbcompliance.pl
#
# This script deals with certain genome browser errors, allowing replacement of 
# the header and gets rid of positions past the reference chromosomes.
#
# Arguments:
#    <see below>
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "usage: gbcompliance.pl
    1) Database (hg19, hg18, rn4)
    2) Input Prefix (Include full path)
    3) Output Prefix (Include full path)
    4) Track Name Prefix
    5) Description Prefix
" unless @ARGV == 5;

my $genome = shift(@ARGV);
my %Chroms;	
if($genome eq "hg18"){
        %Chroms = ("chr1" => '247249719',
                   "chr2" => '242951149',
                   "chr3" => '199501827',
                   "chr4" => '191273063',
                   "chr5" => '180857866',
                   "chr6" => '170899992',
                   "chr7" => '158821424',
                   "chr8" => '146274826',
                   "chr9" => '140273252',
                   "chr10" => '135374737',
                   "chr11" => '134452384',
                   "chr12" => '132349534',
                   "chr13" => '114142980',
                   "chr14" => '106368585',
                   "chr15" => '100338915',
                   "chr16" => '88827254',
                   "chr17" => '78774742',
                   "chr18" => '76117153',
                   "chr19" => '63811651',
                   "chr20" => '62435964',
                   "chr21" => '46944323',
                   "chr22" => '49691432',
                   "chrM" => '16571',
                   "chrX" => '154913754',
                   "chrY" => '57772954',);
}
elsif($genome eq "hg19"){
        %Chroms = ("chr1" => '249250621',
                   "chr2" => '243199373',
                   "chr3" => '198022430',
                   "chr4" => '191154276',
                   "chr5" => '180915260',
                   "chr6" => '171115067',
                   "chr7" => '159138663',
                   "chr8" => '146364022',
                   "chr9" => '141213431',
                   "chr10" => '135534747',
                   "chr11" => '135006516',
                   "chr12" => '133851895',
                   "chr13" => '115169878',
                   "chr14" => '107349540',
                   "chr15" => '102531392',
                   "chr16" => '90354753',
                   "chr17" => '81195210',
                   "chr18" => '78077248',
                   "chr19" => '59128983',
                   "chr20" => '63025520',
                   "chr21" => '48129895',
                   "chr22" => '51304566',
                   "chrM" => '16571',
                   "chrX" => '155270560',
                   "chrY" => '59373566',);
}
elsif($genome eq "rn4"){
        %Chroms = ("chr1" => '267910886',
                   "chr2" => '258207540',
                   "chr3" => '171063335',
                   "chr4" => '187126005',
                   "chr5" => '173096209',
                   "chr6" => '173096209',
                   "chr7" => '143002779',
                   "chr8" => '129041809',
                   "chr9" => '113440463',
                   "chr10" => '110718848',
                   "chr11" => '87759784',
                   "chr12" => '46782294',
                   "chr13" => '111154910',
                   "chr14" => '112194335',
                   "chr15" => '109758846',
                   "chr16" => '90238779',
                   "chr17" => '97296363',
                   "chr18" => '87265094',
                   "chr19" => '59218465',
                   "chrX" => '160699376',);
}
else{die "$genome is not one of the following genomes: {hg18, hg19, rn4}";}

my $inputprefix = shift(@ARGV);	
my $outputprefix = shift(@ARGV);
my $trackname = shift(@ARGV);
my $description = shift(@ARGV);

print "Initializing...\n";

foreach my $chr (sort keys %Chroms)
{
	
	print "Processing Chromosome $chr...\n";
	my $inputfile = $inputprefix . $chr . ".bed";
	open(IN, "<$inputfile") or die "gbcompliance.pl: Error: Missing Chromosome, cannot open $inputfile input BED file";

	my $outputfile = $outputprefix . $chr . ".bed";
	open(OUT, ">$outputfile") or die "gbcompliance.pl: Error: cannot open $outputfile output file";

	print OUT "track name=", $trackname, $chr, " description=", $description, $chr, " useScore=0 itemRgb=On db=", $genome, "\n";

	my $original = <IN>; # Remove header
	my @line = split("\t",$original);
	if ($line[0] =~ /^chr/)
	{
		if ($line[1] <= $Chroms{$chr} && $line[2] <= $Chroms{$chr})
		{
			print OUT $original;
		}
		
	} 
	else {}

	while(<IN>)
	{
		my $original = $_;
		my @line = split("\t",$original);
		if ($line[1] <= $Chroms{$chr} && $line[2] <= $Chroms{$chr})
		{
			print OUT $original;
		}
	}
}




