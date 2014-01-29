#!/usr/bin/perl -w

# Keith Dunaway, 6/8/13 version

# Makes it so you only have to go through the file once.

# Diane Schroeder, 7/15/10 version


# BS Seeker changes the names of the chromosomes to numbers.  For the
#   arguements for this program, use the same BS Seeker number naming 
#   conventions, but feel free to use the actual chromosome name in the
#   file name.

#----------------------------------------------------------------------

# Since the hash uses start position as the key, you can only have one read
#   with a given start position.  For our purposes, this is okay since muliple
#   reads per start position might be clonal (from amplification).  So, we'll
#   keep one (the last) copy of each read at each position in each orientation,
#   letting it overwrite prvious ones.  At the end we'll give statistics about
#   what percent of reads had "clones", average number of clones per position,
#   etc.  For now this won't include orientation, so two reads with the same
#   start site but opposite orientations (and thus not clonal) will overwrite
#   each other.

# For now, this program will only use reads from one chromosome.  The BS Seeker
#   program can have reads mapped to other chromosomes, but they will be
#   disgarded.  Since there will be so much data anyway, this is probably 
#   best, especially if it will be later used for custom tracks in the
#   Santa Cruz genome browser.  Just run the program multiple times, once 
#   for each chromosome.

# Arguments:
#   1) name of BS_Seeker output file
#   2) output file of sorted data
#   3) chromosome (ex. "0006")


unless (@ARGV == 3) {

    print "Arguments:\n\tBS Seeker input file\n\toutput file name for sorted data\n\tgenome (ex. hg18 or mm9)\n\n";
    die;
}
$input = shift(@ARGV);
$outprefix = shift(@ARGV);
$genome = shift(@ARGV);

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
	           '0023' => "chrY",
	           '0023' => "chrM",);
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

my $sortedinput =  $outprefix . "sorted.txt";
my $commandline = "sort -nk4 " . $input . " > " $sortedinput;
`$commandline`;
open (INPUT, "<$sortedinput") or die "Couldn't open $sortedinput\n";

my $chrom = "null";
my $total_reads = 0;
my $total_positions = 0;
my $total_clonal = 0;
my $bshash;

# We'll create a big hash of the data, with the keys being the start postions.
#   The values will be a 2-element array.  The first element will be the 
#   entire BS Seeker line.  The second element will be the count of clones for 
#   that position.  The keys will then be sorted
#   and values written to a file.  Hopefully you have memory....


while($line = <INPUT>){
    chomp $line;
    my @fields = split /\t/, $line;
    $fields[3] =~ m/(\d+)(\+|-)0*(\d+)/;
    my $this_chrom = $1;
    my $this_chrom = $3;
	
	if($linechrom ne $chrom){
		foreach $start ( sort {$a <=> $b} keys %bshash ) {
		    print OUTPUT "$bshash{$start}[0]\n";
		    if ($bshash{$start}[1] > $max_clonal) {
				$max_clonal = $bshash{$start}[1];
		    }
		}
		print "\n\nTotal number of reads on chromosome $chrom:  $total_reads\n\n";
		print "Total number of positions:  $total_positions\n\n";
		my $temp = $total_reads / $total_positions;
		print "Average number of reads per position:  $temp\n\n";
		print "Number of positions that had multiple clones:  $total_clonal\n\n";
		$temp = $total_clonal / $total_positions;
		print "Percent of positions that had multiple clones:  $temp\n\n";
		print "Maximum number of clones per position:  $max_clonal\n\n";
	
		$total_reads = 0;
		$total_positions = 0;
		$total_clonal = 0;
		for (keys %bshash){
        	delete $bshash{$_};
	    }
		close OUT;
		$chrom = $linechrom;
		my $output = $outprefix . $chrom . ".txt";
		open (OUTPUT, ">$output") or die "Couldn't open $output\n";
	}
	
	$total_reads++;
	if (exists($bshash{$position})){
	    if ($bshash{$position}[1] == 1) {
			$total_clonal++;
	    }
	    $bshash{$position}[1]++;
	}
	else {
	    $bshash{$position}[1] = 1; 
	    $total_positions++;
	}
	$bshash{$position}[0] = $line;
	$max_clonal = 1;  #maximum number of clonal reads per position
}
    
    





while () {



    if ($this_chrom eq $chrom) { 


	#The hash keys will be the start position (col4)



    } #end if (read is type we're interested in)

} #end while (going through input file)



# Now write everything out

