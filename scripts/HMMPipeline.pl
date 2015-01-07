#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: HMMPipeline.pl
# Version: 1.4
# Last Updated: 12/28/2014
#
# This is a wrapper that puts has minimal input and automates output. It is NOT
#   memory efficient.
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

# Global paths (need to change for each computer)
my $BED_to_Korf_path = "/data/scratch/programs/Methylation_Domain_Pipeline/Example_Run_Folder/chromFa/BED_to_Korf_input.pl";
my $RunMethylationDomainPipeline_path = "/data/scratch/programs/Methylation_Domain_Pipeline/Example_Run_Folder/RunMethylationDomainPipeline.pl";
my $hg18fastapath = "/data/scratch/genomes/hg18/chrom_fasta/";
my $hg19fastapath = "/data/scratch/genomes/hg19/chrom_fasta/";
my $hg38fastapath = "/data/scratch/genomes/hg38/chrom_fasta/";
my $rn4fastapath = "/data/scratch/genomes/rn4/chrom_fasta/";
my $rn6fastapath = "/data/scratch/genomes/rn6/chrom_fasta/";

# I/O check
die "This script needs the following arguments:
    1) Permeth bed files prefix (leave off chr##.bed)
    2) Experiment Name (prefix for all outfiles)
    3) genome (hg18, hg19, hg38, rn6, or rn4)
    4) HMM model name (ex: Placenta_30_upd.hmm)
" unless @ARGV == 4;

my $bedprefix = shift(@ARGV);
my $experimentname = shift(@ARGV);
my $genome = shift(@ARGV);
my $HMMmodel = shift(@ARGV);

#die "Paramater 5 (zipout) needs to be Y or N, not $zipout" unless ($zipout eq "Y" | $zipout eq "y" | $zipout eq "N" | $zipout eq "n");

my %Chroms;
my $GenomeFASTApath;
my $CPGIslandsbed;
if($genome eq "hg18"){
		$GenomeFASTApath = $hg18fastapath;
		$CPGIslandsbed = "hg18_genome_CGI.bed";
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
                   "chrX" => '154913754',
                   "chrY" => '57772954',);
}
elsif($genome eq "hg19"){
		$GenomeFASTApath = $hg19fastapath;
		$CPGIslandsbed = "hg19_genome_CGI.bed";
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
                   "chrX" => '155270560',
                   "chrY" => '59373566',);
}
elsif($genome eq "hg38"){
		$GenomeFASTApath = $hg38fastapath;
		$CPGIslandsbed = "hg38_genome_CGI.bed";
        %Chroms = ("chr1" => '248956422',
                   "chr2" => '242193529',
                   "chr3" => '198295559',
                   "chr4" => '190214555',
                   "chr5" => '181538259',
                   "chr6" => '170805979',
                   "chr7" => '159345973',
                   "chr8" => '145138636',
                   "chr9" => '138394717',
                   "chr10" => '133797422',
                   "chr11" => '135086622',
                   "chr12" => '133275309',
                   "chr13" => '114364328',
                   "chr14" => '107043718',
                   "chr15" => '101991189',
                   "chr16" => '90338345',
                   "chr17" => '83257441',
                   "chr18" => '80373285',
                   "chr19" => '58617616',
                   "chr20" => '64444167',
                   "chr21" => '46709983',
                   "chr22" => '50818468',
                   "chrM" => '16569',
                   "chrX" => '156040895',
                   "chrY" => '57227415',);
}
elsif($genome eq "rn4"){
		$GenomeFASTApath = $rn4fastapath;
		$CPGIslandsbed = "rn4_genome_CGI.bed";
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
elsif($genome eq "rn6"){
		$GenomeFASTApath = $rn6fastapath;
		$CPGIslandsbed = "rn6_genome_CGI.bed";
        %Chroms = ("chr1" => '282763074',
                   "chr2" => '266435125',
                   "chr3" => '177699992',
                   "chr4" => '184226339',
                   "chr5" => '173707219',
                   "chr6" => '147991367',
                   "chr7" => '145729302',
                   "chr8" => '133307652',
                   "chr9" => '122095297',
                   "chr10" => '112626471',
                   "chr11" => '90463843',
                   "chr12" => '52716770',
                   "chr13" => '114033958',
                   "chr14" => '115493446',
                   "chr15" => '111246239',
                   "chr16" => '90668790',
                   "chr17" => '90843779',
                   "chr18" => '88201929',
                   "chr19" => '62275575',
                   "chr20" => '56205956',
                   "chrX" => '159970021',);
}
else{die "$genome is not hg18, hg19, hg38, rn6, or rn4";}

################
# Main Section #
################

my $commandline = "";
my $HMMcommandline = "perl " . $RunMethylationDomainPipeline_path . " " . $HMMmodel . " " . $CPGIslandsbed;

foreach my $chrom (sort keys %Chroms) {
        my $end = $Chroms{$chrom};
        $commandline = "perl " . $BED_to_Korf_path . " " . $bedprefix . $chrom . ".bed " . $experimentname . "_" . $chrom . ".b2k.fasta " . $GenomeFASTApath . $chrom . ".fa 1 " . $end;
        $HMMcommandline = $HMMcommandline . " " . $experimentname . "_" . $chrom . ".b2k.fasta";
        print $commandline , "\n";
        `$commandline`;
}

# Move new b2k.fasta files to chromFa folder for next step processing
$commandline = "mv " . $experimentname . "* chromFa/.";
print $commandline , "\n";
`$commandline`;

print $HMMcommandline , "\n";
`$HMMcommandline`;

$commandline = "mv HIGH_Final.bed HMD_" . $experimentname . ".bed";
print $commandline , "\n";
`$commandline`;
$commandline = "mv PARTIAL_Final.bed PMD_" . $experimentname . ".bed";
print $commandline , "\n";
`$commandline`;

$commandline = "perl /data/scratch/programs/perl_script/change_singlebedhead.pl PMD_" . $experimentname . ".bed PMD_" . $experimentname;
print $commandline , "\n";
`$commandline`;
$commandline = "perl /data/scratch/programs/perl_script/change_singlebedhead.pl HMD_" . $experimentname . ".bed HMD_" . $experimentname;
print $commandline , "\n";
`$commandline`;


