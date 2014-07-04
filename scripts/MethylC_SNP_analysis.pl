#!/usr/bin/perl 
use strict; use warnings;

##########################################################################################
# Author: Keith Dunaway
# Email: kwdunaway@ucdavis.edu
# Script Name: MethylC_SNP_analysis.pl
# Version: 1.0
# Last Updated: April 1, 2014
#
# This script looks through all raw fastq sequencing reads and finds the reads that
# have the Line1 pattern.  Then, it puts those reads in a separate file.  Also, the
# script quantifies methylation of these sequences across three potential 
# methylation sites.
#
# Arguments:
#   <see below>
# To call:
#   perl /data/scratch/programs/perl_script/Line1_analysisv2.pl Stats_L1v2.txt Sample_JLD018/JLD018_filtered.fq_LINE1reads.fq Sample_JLDS019/JLD019_filtered.fq_LINE1reads.fq Sample_JLKD006/JLKD_006_filtered.fq_LINE1reads.fq Sample_JLD017/JLD017_filtered.fq_LINE1reads.fq Sample_JLKD008/JLKD008_filtered.fq_LINE1reads.fq Sample_JLKD007/JLKD007_filtered.fq_LINE1reads.fq Sample_JLKD002/JLKD002_filtered.fq_LINE1reads.fq Sample_JLKD003/JLKD003_filtered.fq_LINE1reads.fq Sample_JLKD005/JLKD_005_filtered.fq_LINE1reads.fq Sample_JLKD001/JLKD001_filtered.fq_LINE1reads.fq Sample_JLKD004/JLKD_004_filtered.fq_LINE1reads.fq
#
##########################################################################################



####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Sequences to Assay
    2) Maximum identical copies of reads used for analysis (suggest 1)
#    ) Stats outfile (also makes a read file with counts named seqcount_statsoutfile)
#    ) Output suffix (creates an output FASTQ file with only the reads that have LINE1)
    3-?) Input fastq file
" unless @ARGV > 2;

my $seqs_to_assay = shift(@ARGV);	
#Convert inseq to all uppercase so script doesn't break by case issues
$seqs_to_assay =~ tr/[a-z]/[A-Z]/;
my @inseqs = split("," , $seqs_to_assay);
my $maxcopies = shift(@ARGV);
my @infiles = @ARGV;



my @searchterms;
my @captureterms;
my @searchtermswithSNPs;
my @capturetermswithSNPs;
for(my $c = 0; $c < @inseqs; $c++){
	my $inseq = $inseqs[$c];
	if($inseq eq "LINE1"){
		# Line1 sequence used to analyze
		# CTCGTGGTGCGCCGTTTCTTAAGCCGG = original sequence, reverse strand in genome
		# TTYGTGGTGYGTYGTTTTTTAAKTYGG = BS converted sequence
		$inseq = "CTCGTGGTGCGCCGTTTCTTAAGCCGG";

		# Terms with SNP considerations	
#		@searchterms = ("TT[CT][GA]TGGTG[CT][GA]T[CT][GA]TTTTTTAA",
#						"TT[GA][GA][GA][GA][GA][GA][CT][GA][GA][CT][GA]T[GA]TT[GA][CT][GA][GA][GA]");
#		@captureterms = ("TT([CT][GA])TGGTG([CT][GA])T([CT][GA])TTTTTTAAGT([CT][GA])",
#						"([CT][GA])[GA]TTT[GA][GA][GA][GA][GA][GA]([CT][GA])[GA]([CT][GA])T[GA]TT[GA]([CT][GA])[GA][GA]");
	}
	elsif($inseq eq "ALU"){
		$inseq = "TAGCCGGGCGCGGTGGCGGGCG";
		# Terms with SNP consideration
	#	@searchterms = ("TAGT[CT][GA]GG[CT][GA][CT][GA]GTGG[CT][GA]GG[CT][GA]",
	#					"TAAC[CT][GA]AA[CT][GA][CT][GA]ATAA[CT][GA]AA[CT][GA]", 
	#					"[CT][GA]TT[CT][GA]TTAT[CT][GA][CT][GA]TT[CT][GA]GTTA", 
	#					"[CT][GA]CC[CT][GA]CCAC[CT][GA][CT][GA]CC[CT][GA]ACTA");
	}
	else{}
	@searchterms = GetSearchTerms($inseq);
	@captureterms = GetCaptureTerms($inseq);
	@searchtermswithSNPs = GetSearchTerms_withSNPs($inseq);
	@capturetermswithSNPs = GetCaptureTerms_withSNPs($inseq);
}


#print $inseqs[0] , "\n" , $searchterms[0] , "\n" , $searchtermswithSNPs[0] , "\n" , $searchterms[1] , "\n" , $searchtermswithSNPs[1] , "\n\n";


#################
# In Files Loop #
#################

while(@infiles){
	my $infile = shift(@infiles);
#
	print "\nFiltering $infile for all matched reads\n";
	my @matched_reads = GetAllMatchReads($infile, @searchterms);

#	my $outfile = $infile . "_" . $inseqs[0] . "reads.fq";
#	open(OUT, ">$outfile") or die "cannot open $outfile infile";
#	for(my $t = 0; $t < @matched_reads; $t++){
#		print OUT "\n" , $searchterms[$t], "\n\n\n", $matched_reads[$t];
#	}
#	close OUT;
#}





###############
# Subroutines #
###############

sub GetSearchTerms{
    my $dna = shift;
    $dna =~ tr/[a-z]/[A-Z]/;
    my $rcdna = reverse_complement($dna);
    my @results;

	my $bsdna = bs_convert_ss($dna);
	my $bsrcdna = bs_convert_ss($rcdna);

	push(@results, GetSearchTerm($bsdna));
	push(@results, GetSearchTerm($bsrcdna));
	
	return(@results);
}
sub GetSearchTerms_withSNPs{
    my $dna = shift;
    $dna =~ tr/[a-z]/[A-Z]/;
    my $rcdna = reverse_complement($dna);
    my @results;

	my $bsdna = bs_convert_ss_withSNPs($dna);
	my $bsrcdna = bs_convert_ss_withSNPs($rcdna);

	push(@results, GetSearchTerm($bsdna));
	push(@results, GetSearchTerm($bsrcdna));
	
	return(@results);
}
sub GetCaptureTerms{
    my $dna = shift;
    $dna =~ tr/[a-z]/[A-Z]/;
    my $rcdna = reverse_complement($dna);
    my @results;

	my $bsdna = bs_convert_ss($dna);
	my $bsrcdna = bs_convert_ss($rcdna);
	my $search_bsdna = GetSearchTerm($bsdna);
	my $search_bsrcdna = GetSearchTerm($bsrcdna);
	
	push(@results, annotate_cpgs($search_bsdna));
	push(@results, annotate_cpgs($search_bsrcdna));
	
	return(@results);
}
sub GetCaptureTerms_withSNPs{
    my $dna = shift;
    $dna =~ tr/[a-z]/[A-Z]/;
    my $rcdna = reverse_complement($dna);
    my @results;

	my $bsdna = bs_convert_ss_withSNPs($dna);
	my $bsrcdna = bs_convert_ss_withSNPs($rcdna);
	my $search_bsdna = GetSearchTerm($bsdna);
	my $search_bsrcdna = GetSearchTerm($bsrcdna);
	
	push(@results, annotate_cpgs($search_bsdna));
	push(@results, annotate_cpgs($search_bsrcdna));
	
	return(@results);
}
sub reverse_complement {
    my $dna = shift;

	# reverse the DNA sequence
    my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

    return $revcomp;
}
sub bs_convert_ss {
    my $bsdna = shift;

    #Put SNP Y (C or T) for every CpG (because the C's might be methylated
    $bsdna =~ s/CG/YG/g;

    #Convert all non CpGs C's to T's (because they won't be methylated)
    $bsdna =~ s/C/T/g;

    return $bsdna;
}
sub bs_convert_ss_withSNPs {
    my $bsdna = shift;

    #Put SNP Y (C or T) for every CpG (because the C's might be methylated
    $bsdna =~ s/CG/YG/g;

    #Convert all non CpGs C's to T's (because they won't be methylated)
    $bsdna =~ s/C/T/g;

    $bsdna =~ s/A/R/g;
    $bsdna =~ s/G/R/g;

    return $bsdna;
}
sub bs_convert_os {
    my $bsdna = shift;

    #Put SNP R (G or A) for every CpG (because the C's might be methylated on opposite strand
    $bsdna =~ s/CG/CR/g;

    #Convert all non CpGs G's to A's (because the C's won't be methylated)
    $bsdna =~ s/G/A/g;

    return $bsdna;
}
sub bs_convert_os_withSNPs {
    my $bsdna = shift;

    #Put SNP R (G or A) for every CpG (because the C's might be methylated on opposite strand
    $bsdna =~ s/CG/CR/g;

    #Convert all non CpGs G's to A's (because the C's won't be methylated)
    $bsdna =~ s/G/A/g;

    $bsdna =~ s/C/Y/g;
    $bsdna =~ s/T/Y/g;

    return $bsdna;
}
sub GetSearchTerm {
	my $dna = shift;

    my %ambig_codes = (
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
        N => '[ACGT]',
    );
	foreach my $ambig (keys %ambig_codes){
		$dna =~ s/$ambig/$ambig_codes{$ambig}/g;
	}
	return($dna);
}
sub annotate_cpgs{
    my $searchstring = shift;
    $searchstring =~ s/Y/[CT]/g;
    $searchstring =~ s/R/[GA]/g;
    $searchstring =~ s/\[CT\]\[GA\]/\)\([CT][GA]\)\(/g;    
    $searchstring =~ s/\[CT\]G/\([CT][GA]\)\(/g;
    $searchstring = "(" . $searchstring . ")";
    return($searchstring);
}
sub PrintBSseq_both{
	my $dna = shift;
	my $retseq = "Original Sequence:\t" . $dna . "\n";

	my $tmp = bs_convert_ss($dna);
	$retseq = $retseq . "Forward BS conv:  \t" . $tmp . "\n";

	$tmp = reverse_complement($dna);
	$tmp = bs_convert_ss($tmp);
	$retseq = $retseq . "Reverse BS conv:  \t" . $tmp . "\n";

	return($retseq);
}

# Not finished
sub BS_Convert{
	# Takes in a sequence and returns a hash with the following information:
	# 
	#  NOTE: This only gives information for same strand conversion. If
	#        you need opposite strand information, use BS_Convert_All
	
    my $dna = shift;
    
    #Convert everything to uppercase
    $dna =~ tr/[a-z]/[A-Z]/;
	#Get reverse complement
    my $rcdna = reverse_complement($dna);
    
    #Hash with multiple information about sequence to assay in it
    my %sequences;

    #Get bisulfite converted sequence, both strand conversion
    $sequences{"start"}{"seq"} = $dna;
    $sequences{"revcomp"}{"seq"} = $rcdna;
    $sequences{"bsfs"}{"seq"} = bs_convert_ss($dna);
    $sequences{"bsfo"}{"seq"} = bs_convert_os($dna);
    $sequences{"bsrs"}{"seq"} = bs_convert_ss($rcdna);
    $sequences{"bsro"}{"seq"} = bs_convert_os($rcdna);

	#Get search string for CpG's
    $sequences{"bsfs"}{"search"} = annotate_cpgs($sequences{"bsfs"}{"seq"});
    $sequences{"bsfo"}{"search"} = annotate_cpgs($sequences{"bsfo"}{"seq"});
    $sequences{"bsrs"}{"search"} = annotate_cpgs($sequences{"bsrs"}{"seq"});
    $sequences{"bsro"}{"search"} = annotate_cpgs($sequences{"bsro"}{"seq"});
    
    #Find SNP sites
    @{$sequences{"bsfs"}{"SNPs"}} = find_SNPsites($sequences{"bsfs"}{"seq"}, 'Y');
    @{$sequences{"bsfo"}{"SNPs"}} = find_SNPsites($sequences{"bsfo"}{"seq"}, 'R');
    @{$sequences{"bsrs"}{"SNPs"}} = find_SNPsites($sequences{"bsrs"}{"seq"}, 'Y');
    @{$sequences{"bsro"}{"SNPs"}} = find_SNPsites($sequences{"bsro"}{"seq"}, 'R');

#    my %sequences = ($dna_bsfs,$dna_bsfo,$dna_bsrs,$dna_bsro);
	return(%sequences);
}

#Returns SNP locations as an array of numbers
sub find_SNPsites{
	# String containing sequence to be looked through for a SNP
	my $string = shift;
	# Character to be looked through the string for, most likely Y for R.
	my $char = shift;
	# Resulting locations to be returned
	my @locations;

	my $offset = 0;
	my $result = index($string, $char, $offset);
	while ($result != -1) {
		push(@locations, $result);
		$offset = $result + 1;
    	$result = index($string, $char, $offset);
  }
	
	return(@locations);
}


sub Print_BS_Sequences {
	my ($seq_ref) = @_;
#    my $params = shift;
    my %sequences = %$seq_ref;
    
	my $print_statement =
	"Original:             \t" . $sequences{"start"}{"seq"} . "\n" .
	"BS For Same:          \t" . $sequences{"bsfs"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsfs"}{"search"} . "\n\n" .
	
	"Original:             \t" . $sequences{"start"}{"seq"} . "\n" .
	"BS For Opp:           \t" . $sequences{"bsfo"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsfo"}{"search"} . "\n\n" .

	"Reverse:              \t" . $sequences{"revcomp"}{"seq"} . "\n" .
	"BS Rev Same:          \t" . $sequences{"bsrs"}{"seq"} . "\n" .
	"Searchstring:         \t" . $sequences{"bsrs"}{"search"} . "\n\n" .
	
	"Reverse:              \t" . $sequences{"revcomp"}{"seq"} . "\n" .
	"BS Rev Opp:           \t" . $sequences{"bsro"}{"seq"} . "\n" . 
	"Searchstring:         \t" . $sequences{"bsro"}{"search"} . "\n\n";
	
	return($print_statement);
}

sub GetAllSNP {
    my $infile = shift;
    my $maxcopies = shift;
	open(IN, "<$infile") or die "cannot open $infile infile";
	my ($seq_ref) = @_;
    my %sequences = %$seq_ref;
    my %SNPhash;
    
    my $counter = 0;

	while (<IN>) {
		$counter++;
		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;	    

		foreach my $key (keys %sequences) {
	    	if($seq =~ /$sequences{$key}{"search"}/){
				my $match = $1;
				my $searchterm = $sequences{$key}{"search"};
				if (defined $SNPhash{$seq}){
					$SNPhash{$key}{$seq}{"copy_number"}++;
				}
				else{
					$SNPhash{$key}{$seq}{"copy_number"} = 1;
					@{$SNPhash{$key}{$seq}{"SNPs"}} = FindSeqSNPs($searchterm, @{$sequences{$key}{"SNPs"}});
				}
			}
		}
	}
	close IN;
	return(%SNPhash);
}

sub GetAllMatchReads {
    my $infile = shift;
	open(IN, "<$infile") or die "cannot open $infile infile";
	my @sequences_array = @_;
    my @matched_reads;
    my @counts;
    for (my $r = 0; $r < @sequences_array; $r++){
		$matched_reads[$r] = "";
		$counts[$r] = 0;
    }
    my $counter = 0;
	while (<IN>) {
		$counter++;
		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
	    		$matched_reads[$n] = $matched_reads[$n] . $ID . $seq . "\n" . $third . $quality;
	    		$counts[$n]++;
			}
		}
	}
	close IN;
	
	for(my $n = 0; $n < @sequences_array; $n++){ print "$sequences_array[$n]\t\t$counts[$n]\n"; }
	
	return(@matched_reads);
}

sub CaptureSNPs {
    my $infile = shift;
	open(IN, "<$infile") or die "cannot open $infile infile";
	my @sequences_array = @_;
    my %SNPsequences;
#    my @counts;
#    for (my $r = 0; $r < @sequences_array; $r++){
#		$matched_reads[$r] = "";
#		$counts[$r] = 0;
#    }
#    my $counter = 0;
	while (<IN>) {
#		$counter++;
#		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
	    		if(defined $SNPsequences{$sequences_array[$n]}{$1}){ $SNPsequences{$sequences_array[$n]}{$1}++; }
	    		else{ $SNPsequences{$sequences_array[$n]}{$1} = 1; }
#	    		$matched_reads[$n] = $matched_reads[$n] . $ID . $seq . "\n" . $third . $quality;
#	    		$counts[$n]++;
			}
		}
	}
	close IN;
	
#	for(my $n = 0; $n < @sequences_array; $n++){ print "$sequences_array[$n]\t\t$counts[$n]\n"; }
	
	return(%SNPsequences);
}

sub CaptureL1SNPs {
    my $infile = shift;
	open(IN, "<$infile") or die "cannot open $infile infile";
	my @sequences_array = @_;
    my %SNPsequences;
    $SNPsequences{$sequences_array[0]}{"count"} = 0;
    $SNPsequences{$sequences_array[0]}{1}{"CG"} = 0;
    $SNPsequences{$sequences_array[0]}{1}{"TG"} = 0;
    $SNPsequences{$sequences_array[0]}{1}{"CA"} = 0;
    $SNPsequences{$sequences_array[0]}{1}{"TA"} = 0;
    $SNPsequences{$sequences_array[0]}{2}{"CG"} = 0;
    $SNPsequences{$sequences_array[0]}{2}{"TG"} = 0;
    $SNPsequences{$sequences_array[0]}{2}{"CA"} = 0;
    $SNPsequences{$sequences_array[0]}{2}{"TA"} = 0;
    $SNPsequences{$sequences_array[0]}{3}{"CG"} = 0;
    $SNPsequences{$sequences_array[0]}{3}{"TG"} = 0;
    $SNPsequences{$sequences_array[0]}{3}{"CA"} = 0;
    $SNPsequences{$sequences_array[0]}{3}{"TA"} = 0;
    $SNPsequences{$sequences_array[0]}{4}{"CG"} = 0;
    $SNPsequences{$sequences_array[0]}{4}{"TG"} = 0;
    $SNPsequences{$sequences_array[0]}{4}{"CA"} = 0;
    $SNPsequences{$sequences_array[0]}{4}{"TA"} = 0;
    $SNPsequences{$sequences_array[1]}{"count"} = 0;
    $SNPsequences{$sequences_array[1]}{1}{"CG"} = 0;
    $SNPsequences{$sequences_array[1]}{1}{"TG"} = 0;
    $SNPsequences{$sequences_array[1]}{1}{"CA"} = 0;
    $SNPsequences{$sequences_array[1]}{1}{"TA"} = 0;
    $SNPsequences{$sequences_array[1]}{2}{"CG"} = 0;
    $SNPsequences{$sequences_array[1]}{2}{"TG"} = 0;
    $SNPsequences{$sequences_array[1]}{2}{"CA"} = 0;
    $SNPsequences{$sequences_array[1]}{2}{"TA"} = 0;
    $SNPsequences{$sequences_array[1]}{3}{"CG"} = 0;
    $SNPsequences{$sequences_array[1]}{3}{"TG"} = 0;
    $SNPsequences{$sequences_array[1]}{3}{"CA"} = 0;
    $SNPsequences{$sequences_array[1]}{3}{"TA"} = 0;
    $SNPsequences{$sequences_array[1]}{4}{"CG"} = 0;
    $SNPsequences{$sequences_array[1]}{4}{"TG"} = 0;
    $SNPsequences{$sequences_array[1]}{4}{"CA"} = 0;
    $SNPsequences{$sequences_array[1]}{4}{"TA"} = 0;
    
#    my @counts;
#    for (my $r = 0; $r < @sequences_array; $r++){
#		$matched_reads[$r] = "";
#		$counts[$r] = 0;
#    }
#    my $counter = 0;
	while (<IN>) {
#		$counter++;
#		if($counter % 100000 == 0) {print $counter , "\n";}
	    my $ID = $_;
   		my $seq = <IN>;
	  	chop($seq); #gets rid of return character at end of sequence
	    my $third = <IN>;
	    my $quality = <IN>;
		for(my $n = 0; $n < @sequences_array; $n++){
	    	if($seq =~ /$sequences_array[$n]/){
				$SNPsequences{$sequences_array[$n]}{"count"}++;
				$SNPsequences{$sequences_array[$n]}{1}{$1} = $SNPsequences{$sequences_array[$n]}{1}{$1} + 1;
				$SNPsequences{$sequences_array[$n]}{2}{$1} = $SNPsequences{$sequences_array[$n]}{2}{$2} + 1;
				$SNPsequences{$sequences_array[$n]}{3}{$1} = $SNPsequences{$sequences_array[$n]}{3}{$3} + 1;
				$SNPsequences{$sequences_array[$n]}{4}{$1} = $SNPsequences{$sequences_array[$n]}{4}{$1} + 1;
			}
		}
	}
	close IN;
	return(%SNPsequences);
}

sub PrintL1SNPs{
	my ($seq_ref) = @_;
    my %SNPsequences = %$seq_ref;
    my $printstring = "";
	foreach my $searchstring (keys %SNPsequences) {
		$printstring = $printstring . "Reads matched to " . $searchstring . ":\t" . $SNPsequences{$searchstring}{"count"} . "\n";
		for(my $SNPloc = 1; $SNPloc < 5; $SNPloc++){
			$printstring = $printstring . "For CpG " . $SNPloc . ":" . 
			"\tCG:" . $SNPsequences{$searchstring}{$SNPloc}{"CG"} . 
			" (" .  sprintf("%.2f", 100*$SNPsequences{$searchstring}{$SNPloc}{"CG"}/$SNPsequences{$searchstring}{"count"}) . "%)" .
			"\tTG:" . $SNPsequences{$searchstring}{$SNPloc}{"TG"} . 
			" (" .  sprintf("%.2f", 100*$SNPsequences{$searchstring}{$SNPloc}{"TG"}/$SNPsequences{$searchstring}{"count"}) . "%)" .
			"\t\tCA:" . $SNPsequences{$searchstring}{$SNPloc}{"CA"} . 
			" (" .  sprintf("%.2f", 100*$SNPsequences{$searchstring}{$SNPloc}{"CA"}/$SNPsequences{$searchstring}{"count"}) . "%)" .
			"\tTA:" . $SNPsequences{$searchstring}{$SNPloc}{"TA"} .
			" (" .  sprintf("%.2f", 100*$SNPsequences{$searchstring}{$SNPloc}{"TA"}/$SNPsequences{$searchstring}{"count"}) . "%)\n";
		}
		$printstring = $printstring . "\n";
	}
	return($printstring);
}

# returns a hash of characters at specific location in a string
sub FindSeqSNPs {
	my $string = shift;
	my @SNPsites = @_;
	my @Results;
	while(@SNPsites){
		my $char = shift(@SNPsites);
		push (@Results ,substr($string, $char, 1));
	}
	return(@Results);
}
	
sub ConsolidateBSSNPs{
	my ($seq_ref) = @_;
    my %SNPhash = %$seq_ref;
    #my @names = ("bsfs","bsfo","bsrs","bsro");
    my %running_results = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsfs"}});
    my %bsfs = %running_results;
    
#    print "bsfs           \tT=" , $bsfs{"0"}{"T"}, " C=" , $bsfs{"0"}{"C"} , "\tT=", $bsfs{"1"}{"T"}, " C=" , $bsfs{"1"}{"C"} , "\tT=", $bsfs{"2"}{"T"}, " C=" , $bsfs{"2"}{"C"} , "\tT=", $bsfs{"3"}{"T"}, " C=" , $bsfs{"3"}{"C"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";
    
    #BS conversion in forward orientation on opposite strand
    my %bsfo = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsfo"}});
	my $n = -1;
    while(1){
    	$n++;
    	if(defined $running_results{$n}){
	   		if(defined $bsfo{$n}{"C"} or $bsfo{$n}{"T"}){
				foreach my $key (keys %{$bsfo{$n}}) {
				    if(defined $running_results{$n}{$key}){
		    			$running_results{$n}{$key} = $running_results{$n}{$key} + $bsfo{$n}{$key};
	    			}
	   				else {$running_results{$n}{$key} = $bsfo{$n}{$key};}
				}
			}
    		else{
		    	if(defined $running_results{$n}{"C"} && defined $bsfo{$n}{"G"}){
	    			$running_results{$n}{"C"} = $running_results{$n}{"C"} + $bsfo{$n}{"G"};
	    		}
	    		else {$running_results{$n}{"C"} = $bsfo{$n}{"G"};}
		    	if(defined $running_results{$n}{"T"} && defined $bsfo{$n}{"A"}){
	   	 			$running_results{$n}{"T"} = $running_results{$n}{"T"} + $bsfo{$n}{"A"};
	    		}
	    		else {$running_results{$n}{"T"} = $bsfo{$n}{"A"};}
	    	}
	    }
    	else{last;} 
    }
    my $SNPnum = $n;
#    print "bsfo           \tT=" , $bsfo{"0"}{"A"}, " C=" , $bsfo{"0"}{"G"} , "\tT=", $bsfo{"1"}{"A"}, " C=" , $bsfo{"1"}{"G"} , "\tT=", $bsfo{"2"}{"A"}, " C=" , $bsfo{"2"}{"G"} , "\tT=", $bsfo{"3"}{"A"}, " C=" , $bsfo{"3"}{"G"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";

    #BS conversion in reverse orientation on same strand 
    my %bsrs = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsrs"}});
	$n = -1;
	my $t = $SNPnum;
    while(1){
    	$n++;
    	$t--;
    	if(defined $running_results{$t}){
    		if(defined $bsrs{$n}{"A"} or $bsrs{$n}{"G"}){
				foreach my $key (keys %{$bsrs{$n}}) {
			    	if(defined $running_results{$t}{$key}){
		    			$running_results{$t}{$key} = $running_results{$t}{$key} + $bsrs{$n}{$key};
	    			}
	    			else {$running_results{$t}{$key} = $bsrs{$n}{$key};}
				}
   			}
    		else{
		    	if(defined $running_results{$t}{"C"} && defined $bsrs{$n}{"C"}){
	    			$running_results{$t}{"C"} = $running_results{$t}{"C"} + $bsrs{$n}{"C"};
	    		}
	    		else {$running_results{$t}{"C"} = $bsrs{$n}{"C"};}
		    	if(defined $running_results{$t}{"T"} && defined $bsrs{$n}{"T"}){
		    		print $running_results{$t}{"T"} , "\n";
		    		print $bsrs{$n}{"T"} , "\n";
	   	 			$running_results{$t}{"T"} = $running_results{$t}{"T"} + $bsrs{$n}{"T"};
	    		}
	    		else {$running_results{$t}{"T"} = $bsrs{$n}{"T"};}
	    	}
	    }
    	else{last;} 
    }
#    print "bsrs           \tT=" , $bsrs{"3"}{"T"}, " C=" , $bsrs{"3"}{"C"} , "\tT=", $bsrs{"2"}{"T"}, " C=" , $bsrs{"2"}{"C"} , "\tT=", $bsrs{"1"}{"T"}, " C=" , $bsrs{"1"}{"C"} , "\tT=", $bsrs{"0"}{"T"}, " C=" , $bsrs{"0"}{"C"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";
    
    #BS conversion in reverse orientation on opposite strand 
    my %bsro = ConsolidateSNPs($maxcopies , \%{$SNPhash{"bsro"}});
	$n = -1;
	$t = $SNPnum;
    while(1){
    	$n++;
    	$t--;
    	if(defined $running_results{$t}){
    		if(defined $bsro{$n}{"C"} or $bsro{$n}{"T"}){
				foreach my $key (keys %{$bsro{$n}}) {
			    	if(defined $running_results{$t}{$key}){
		    			$running_results{$t}{$key} = $running_results{$t}{$key} + $bsro{$n}{$key};
	    			}
	    			else {$running_results{$t}{$key} = $bsro{$n}{$key};}
				}
   		}
    		else{
		    	if(defined $running_results{$t}{"C"} && defined $bsro{$n}{"G"}){
	    			$running_results{$t}{"C"} = $running_results{$t}{"C"} + $bsro{$n}{"G"};
	    		}
	    		else {$running_results{$t}{"C"} = $bsro{$n}{"G"};}
		    	if(defined $running_results{$t}{"T"} && defined $bsro{$n}{"A"}){
	   	 			$running_results{$t}{"T"} = $running_results{$t}{"T"} + $bsro{$n}{"A"};
	    		}
	    		else {$running_results{$t}{"T"} = $bsro{$n}{"A"};}
	    	}
	    }
   	else{last;} 
    }
#    print "bsfo           \tT=" , $bsrs{"3"}{"A"}, " C=" , $bsrs{"3"}{"G"} , "\tT=", $bsrs{"2"}{"A"}, " C=" , $bsrs{"2"}{"G"} , "\tT=", $bsrs{"1"}{"A"}, " C=" , $bsrs{"1"}{"G"} , "\tT=", $bsrs{"0"}{"A"}, " C=" , $bsrs{"0"}{"G"} , "\n";
#    print "running_results\tT=" , $running_results{"0"}{"T"}, " C=" , $running_results{"0"}{"C"} , "\tT=", $running_results{"1"}{"T"}, " C=" , $running_results{"1"}{"C"} , "\tT=", $running_results{"2"}{"T"}, " C=" , $running_results{"2"}{"C"} , "\tT=", $running_results{"3"}{"T"}, " C=" , $running_results{"3"}{"C"} , "\n";

    my %results;
    %{$results{"total"}} = %running_results;
    %{$results{"bsfs"}} = %bsfs;
    %{$results{"bsfo"}} = %bsfo;
    %{$results{"bsrs"}} = %bsrs;
    %{$results{"bsro"}} = %bsro;
	return(%results);
}

# Returns a hash with the format SNPresults{numberSNPposition}{SNPletter}=count
sub ConsolidateSNPs{
	#   Example of format:
	#   SNPresults{0}{C} = 12
	#   SNPresults{0}{T} = 108
	# This would be 10% C and 90% T

	my $maxcopies = shift;
	my ($seq_ref) = @_;
    my %Readhash = %$seq_ref;
    my %SNPresults;

	foreach my $key (keys %Readhash) {

#		print $Readhash{$key}{"copy_number"} , "\n";

#		my $copiesused = $Readhash{$key}{"copy_number"};
#		if($copiesused > $maxcopies){
			my $copiesused = $maxcopies;
#		}
		my @tmparray = @{$Readhash{$key}{"SNPs"}};
		for(my $n = 0; $n < @tmparray; $n++){
			if(defined $SNPresults{$n}{$tmparray[$n]}){
				$SNPresults{$n}{$tmparray[$n]} = $SNPresults{$n}{$tmparray[$n]} + $copiesused;		
			}
			else{
				$SNPresults{$n}{$tmparray[$n]} = $copiesused;					
			}
		}
	}
	return(%SNPresults);
}


sub PrintSNPs{
	my ($seq_ref) = @_;
    my %results = %$seq_ref;
    
    my $printline = "";

#    my %printresults;    
#	foreach my $file (keys %results) {
#		foreach my $type (keys %{results{$file}}) {
#			%{$printresults{$type}{$file}} = %{$results{$file}{$type}}
#		}
#   }

	foreach my $file (keys %results) {
		foreach my $type (keys %{$results{$file}}) {
			$printline  = $printline . $file . "\t" . $type . "\t";
			foreach my $SNP (keys %{$results{$file}{$type}}) {
			   $printline = $printline . $SNP . "=" . $results{$file}{$type}{$SNP} . " ";
			 }
			 $printline = $printline . "\n";
		}
		$printline = $printline . "\n";
   	}    
	return($printline);
}

sub BS_Convert_All{
    my $dna = shift;
    
    #Convert everything to uppercase
    $dna =~ tr/[a-z]/[A-Z]/;
	#Get reverse complement
    my $rcdna = reverse_complement($dna);
    
    #Hash with multiple information about sequence to assay in it
    my %sequences;

    #Get bisulfite converted sequence, both strand conversion
    $sequences{"start"}{"seq"} = $dna;
    $sequences{"revcomp"}{"seq"} = $rcdna;
    $sequences{"bsfs"}{"seq"} = bs_convert_ss($dna);
    $sequences{"bsfo"}{"seq"} = bs_convert_os($dna);
    $sequences{"bsrs"}{"seq"} = bs_convert_ss($rcdna);
    $sequences{"bsro"}{"seq"} = bs_convert_os($rcdna);

	#Get search string for CpG's
    $sequences{"bsfs"}{"search"} = annotate_cpgs($sequences{"bsfs"}{"seq"});
    $sequences{"bsfo"}{"search"} = annotate_cpgs($sequences{"bsfo"}{"seq"});
    $sequences{"bsrs"}{"search"} = annotate_cpgs($sequences{"bsrs"}{"seq"});
    $sequences{"bsro"}{"search"} = annotate_cpgs($sequences{"bsro"}{"seq"});
    
    #Find SNP sites
    @{$sequences{"bsfs"}{"SNPs"}} = find_SNPsites($sequences{"bsfs"}{"seq"}, 'Y');
    @{$sequences{"bsfo"}{"SNPs"}} = find_SNPsites($sequences{"bsfo"}{"seq"}, 'R');
    @{$sequences{"bsrs"}{"SNPs"}} = find_SNPsites($sequences{"bsrs"}{"seq"}, 'Y');
    @{$sequences{"bsro"}{"SNPs"}} = find_SNPsites($sequences{"bsro"}{"seq"}, 'R');

#    my %sequences = ($dna_bsfs,$dna_bsfo,$dna_bsrs,$dna_bsro);
	return(%sequences);
}

__END__

####################################################################
# Command Line Error Checking. Global Variables and I/O Initiation #
####################################################################

die "This script needs the following arguments:
    1) Sequence to Assay

#    ) Stats outfile (also makes a read file with counts named seqcount_statsoutfile)
#    ) Output suffix (creates an output FASTQ file with only the reads that have LINE1)
    2-?) Input fastq file
" unless @ARGV > 1;

my $inseq = shift(@ARGV);	
#Convert inseq to all uppercase so script doesn't break by case issues
$inseq =~ tr/[a-z]/[A-Z]/;
my @infiles = @ARGV;

# Hard coded for now, but should be able to be set as parameter
my $maxcopies = 1;

#my $stats_filename = shift(@ARGV);	
#my $out_suffix = shift(@ARGV);

#open(STATS, ">$stats_filename") or die "cannot open $stats_filename outfile";
#my $seqcount = "seqcount_" . $stats_filename;
#open(SEQOUT, ">$seqcount") or die "cannot open $seqcount outfile";

# Line1 sequence used to analyze
# TTYGTGGTGYGTYGTTTTTTAAKTYGGTT
# ctcgtggtgcgccgtttcttaagccggtctg
# CTCGTGGTGCGCCGTTTCTTAAGCCGGTC
# CTCGTGGTGCGCCGTTTCTTAA

my @searchterms;
my @captureterms;
if($inseq eq "LINE1"){
#	$inseq = "CTCGTGGTGCGCCGTTTCTTAAGCCGGTC";
#	$inseq = "CTCGTGGTGCGCCGTTTCTTAA";
	$inseq = "CTCGTGGTGCGCCGTTTCTTAAGCCG";
	# Terms without SNP consideration
#	@searchterms = ("TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA",
#					"CTC[GA]TAATAC[GA]CC[GA]TTTCTTAA", 
#					"TTAAGAAA[CT]GG[CT]GTATTA[CT]GAG", 
#					"TTAAAAAAC[GA]AC[GA]CACCAC[GA]AA");

	# Terms with SNP consideration
	
#	@searchterms = ("TT[CT][GA]TGGTG[CT][GA]T[CT][GA]TTTTTTAA",
#					"CT[CT][GA]TAATA[CT][GA]C[CT][GA]TTTCTTAA", 
#					"TTAAGAAA[CT][GA]G[CT][GA]TATTA[CT][GA]AG", 
#					"TTAAAAAA[CT][GA]A[CT][GA]CACCA[CT][GA]AA");
#	@captureterms = ("TT([CT][GA])TGGTG([CT][GA])T([CT][GA])TTTTTTAA([GT])T([CT][GA])",
#					 "CT([CT][GA])TAATA([CT][GA])C([CT][GA])TTTCTTAA([AT])C([CT][GA])", 
#					 "([CT][GA])G([TA])TTAAGAAA([CT][GA])G([CT][GA])TATTA([CT][GA])AG", 
#					 "([CT][GA])A([CA])TTAAAAAA([CT][GA])A([CT][GA])CACCA([CT][GA])AA");
	@searchterms = ("TT[CT][GA]TGGTG[CT][GA]T[CT][GA]TTTTTTAA",
					"CT[CT][GA]TAATA[CT][GA]C[CT][GA]TTTCTTAA", 
					"TT[GA][GA][GA][GA][GA][GA][CT][GA][GA][CT][GA]T[GA]TT[GA][CT][GA][GA][GA]",
					"TTAAAAAA[CT][GA]A[CT][GA]CACCA[CT][GA]AA");
	@captureterms = ("TT([CT][GA])TGGTG([CT][GA])T([CT][GA])TTTTTTAAGT([CT][GA])",
					"([CT][GA])[GA]TTT[GA][GA][GA][GA][GA][GA]([CT][GA])[GA]([CT][GA])T[GA]TT[GA]([CT][GA])[GA][GA]");

}
elsif($inseq eq "ALU"){
	$inseq = "TAGCCGGGCGCGGTGGCGGGCG";
	# Terms with SNP consideration
	@searchterms = ("TAGT[CT][GA]GG[CT][GA][CT][GA]GTGG[CT][GA]GG[CT][GA]",
					"TAAC[CT][GA]AA[CT][GA][CT][GA]ATAA[CT][GA]AA[CT][GA]", 
					"[CT][GA]TT[CT][GA]TTAT[CT][GA][CT][GA]TT[CT][GA]GTTA", 
					"[CT][GA]CC[CT][GA]CCAC[CT][GA][CT][GA]CC[CT][GA]ACTA");
}

my %bs_seq = BS_Convert_All($inseq);

my $printseq = Print_BS_Sequences(\%bs_seq);
print $printseq;


#################
# In Files Loop #
#################

my %results;
my @names = ("bsfs","bsfo","bsrs","bsro");

#my %search_reads;
#for(my $t = 0; $t < @names; $t++){
#	%{$search_reads{$names[$t]}} = %{$bs_seq{$names[$t]}};
#}

#for(my $t = 0; $t < @names; $t++){
#	push(@searchterms, $bs_seq{$names[$t]}{"search"});
#}


#Original:             	CTCGTGGTGCGCCGTTTCTTAA
#BS For Same:          	TT Y  GTGGTG Y  GT Y  GTTTTTTAA
#Searchstring:         	TT[CT]GTGGTG[CT]GT[CT]GTTTTTTAA

#BS For Opp:           	CTC R  TAATAC R  CC R  TTTCTTAA
#Searchstring:         	CTC[GA]TAATAC[GA]CC[GA]TTTCTTAA

#Reverse:              	TTAAGAAACGGCGCACCACGAG
#BS Rev Same:          	TTAAGAAA Y  GG Y  GTATTA Y  GAG
#Searchstring:         	TTAAGAAA[CT]GG[CT]GTATTA[CT]GAG

#BS Rev Opp:           	TTAAAAAAC R  AC R  CACCAC R  AA
#Searchstring:         	TTAAAAAAC[GA]AC[GA]CACCAC[GA]AA


while(@infiles){
	my $infile = shift(@infiles);

	print "\nFiltering $infile for all matched reads\n";
#	my @matched_reads = GetAllMatchReads($infile, @searchterms);
#	my $outfile = $infile . "_l1reads.fq ";
#	open(OUT, ">$outfile") or die "cannot open $outfile infile";
#	for(my $t = 0; $t < @matched_reads; $t++){
#		print OUT "\n" , $searchterms[$t], "\n\n\n", $matched_reads[$t];
#	}
	close OUT;

#	my %SNPsHash = CaptureSNPs($outfile, @captureterms);

	my %SNPsHash = CaptureL1SNPs($infile, @captureterms);
	my $pstring = PrintL1SNPs(\%SNPsHash);
	print $pstring;

#	my %SNPsHash = CaptureSNPs($infile, @captureterms);
#	my $cutoff = 2;
#	foreach my $searchterm (keys %SNPsHash){
#		my $total = 0;
#		my $first = "0";
#		foreach my $key (sort { $SNPsHash{$searchterm}{$b} <=> $SNPsHash{$searchterm}{$a} } keys %{$SNPsHash{$searchterm}}){
#			if($first eq "0") { 
#				$first = $key; 
#				print $key , "\t" , $SNPsHash{$searchterm}{$key}, "\n";
#			}
#			elsif($SNPsHash{$searchterm}{$key} >= $cutoff){
#				my $s1 = $first;
#				my $s2 = $key;	
#				my @s1 = split //,$s1;
#				my @s2 = split //,$s2;
#				my $i = 0;
#				foreach  (@s1) {
#				    if ($_ ne $s2[$i]) {print "$s2[$i]";}
#				    else {print " ";}
#			    	$i++;
#				}
#	       		print "\t" , $SNPsHash{$searchterm}{$key}, "\n";
#		    }
#				$total += $SNPsHash{$searchterm}{$key};


#			if ($SNPsHash{$searchterm}{$key} > 1){
#				print $key , "\t" , $SNPsHash{$searchterm}{$key} , "\n";
#				$total += $SNPsHash{$searchterm}{$key};
#			}
#		}
#		print "Total:\t$total\n\n";
#	}



#	print "\nStarting file: $infile\n";
#	%{$results{$infile}} = GetAllSNP($infile, $maxcopies, \%search_reads);
#	print "Consolidating BSSNPs results\n";
#	%{$results{"processed"}{$infile}} = ConsolidateBSSNPs(\%{$results{$infile}});
}

#$printseq = PrintSNPs(\%{$results{"processed"}});
#print $printseq;




