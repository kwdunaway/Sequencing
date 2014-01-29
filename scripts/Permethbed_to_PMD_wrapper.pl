#!/usr/bin/perl
use warnings;
use strict;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Last Edit Date: 12-5-2013
# Script Name: Permethbed_to_PMD_wrapper.pl
#
# Wrapper that uses other scripts to take a folder of permeth bed files and call
# PMDs/HMDs using a designated model.
#
# Requires the following scripts:
#     BED_to_Korf_input.pl
# 
# Arguments:
#    1) Output prefix
#    2) Input CpG File GTF (cpgIslandExt.gtf)
#    3-?) Input fasta file from perl script BED_to_Korf_input.pl 
#
##########################################################################################


# location of fasta files for different genomes

$hg18fa = "";
my %Chroms;
if($genome eq "hg18"){
	%Chroms = ('0001' => "chr1",
	           '0002' => "chr2",
	           '0003' => "chr3",
	           '0004' => "chr4",
	           '0005' => "chr5",
	           '0006' => "chr6",
	           '0007' => "chr7",
	           '0008' => "chr8",
	           '0009' => "chr9",
	           '0010' => "chr10",
	           '0011' => "chr11",
	           '0012' => "chr12",
	           '0013' => "chr13",
	           '0014' => "chr14",
	           '0015' => "chr15",
	           '0016' => "chr16",
	           '0017' => "chr17",
	           '0018' => "chr18",
	           '0019' => "chr19",
	           '0020' => "chr20",
	           '0021' => "chr21",
	           '0022' => "chr22",
	           '0023' => "chrX",
	           '0024' => "chrY",
	           '0025' => "chrM",);
}
elsif($genome eq "mm9"){
	%Chroms = ('0011' => "chr1",
	           '0012' => "chr2",
	           '0013' => "chr3",
	           '0014' => "chr4",
	           '0015' => "chr5",
	           '0016' => "chr6",
	           '0017' => "chr7",
	           '0018' => "chr8",
	           '0019' => "chr9",
	           '0001' => "chr10",
	           '0002' => "chr11",
	           '0003' => "chr12",
	           '0004' => "chr13",
	           '0005' => "chr14",
	           '0006' => "chr15",
	           '0007' => "chr16",
	           '0008' => "chr17",
	           '0009' => "chr18",
	           '0010' => "chr19",
	           '0021' => "chrX",
	           '0022' => "chrY",
	           '0020' => "chrM",);
}
else{die "$genome is not hg18 or mm9";}

__END__
hg18:
chr1	1	247249719
chr2	1	242951149
chr3	1	199501827
chr4	1	191273063
chr5	1	180857866
chr6	1	170899992
chr7	1	158821424
chr8	1	146274826
chr9	1	140273252
chr10	1	135374737
chr11	1	134452384
chr12	1	132349534
chr13	1	114142980
chr14	1	106368585
chr15	1	100338915
chr16	1	88827254
chr17	1	78774742
chr18	1	76117153
chr19	1	63811651
chr20	1	62435964
chr21	1	46944323
chr22	1	49691432
chrX	1	154913754
chrY	1	57772954


