#!/usr/bin/perl -w

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


unless ($ARGV[2]) {

    print "Arguments:\n\tBS Seeker input file\n\toutput file name for sorted data\n\tchromosome (ex. 0006)\n\n";
    die;
}

open (INPUT, "<$ARGV[0]")
    or die "Couldn't open $ARGV[0]\n";

open (OUTPUT, ">$ARGV[1]")
    or die "Couldn't open $ARGV[2]\n";

$chrom = $ARGV[2];

$total_reads = 0;
$total_positions = 0;
$total_clonal = 0;

# We'll create a big hash of the data, with the keys being the start postions.
#   The values will be a 2-element array.  The first element will be the 
#   entire BS Seeker line.  The second element will be the count of clones for 
#   that position.  The keys will then be sorted
#   and values written to a file.  Hopefully you have memory....


while ($line = <INPUT>) {
    chomp $line;
    @fields = split /\t/, $line;
	# everything before the + or - is the chromosome name
    $fields[3] =~ m/([a-zA-Z0-9_]*)(\+|-)0*(\d+)/;
    $this_chrom = $1;
    $position = $3;
#	die "$fields[3] \n$this_chrom\n$chrom\n$position\n";
    if ($this_chrom eq $chrom) { 
		$total_reads++;
		#The hash keys will be the start position (col4)
		if (exists($bshash{$position})) {
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

    } #end if (read is type we're interested in)

} #end while (going through input file)


$max_clonal = 1;  #maximum number of clonal reads per position

# Now write everything out

foreach $start ( sort {$a <=> $b} keys %bshash ) {

    print OUTPUT "$bshash{$start}[0]\n";

    if ($bshash{$start}[1] > $max_clonal) {

	$max_clonal = $bshash{$start}[1];
    }
}

print "\n\nTotal number of reads on chromosome $chrom:  $total_reads\n\n";
print "Total number of positions:  $total_positions\n\n";

$temp = $total_reads / $total_positions;

print "Average number of reads per position:  $temp\n\n";

print "Number of positions that had multiple clones:  $total_clonal\n\n";

$temp = $total_clonal / $total_positions;

print "Percent of positions that had multiple clones:  $temp\n\n";

print "Maximum number of clones per position:  $max_clonal\n\n";
