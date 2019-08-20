# Combover

Create a "comb" of which positions to record split kmers in a read set.  You choose how far apart the tines are with `-k`.

This repo is just a proof of principle only.

## Usage

First, "sketch" any set of fastq files, one at a time.

    Combover: sketches a fastq file for split kmers and variant
      scripts/comboverSketch.pl -k 12 file.fastq.gz > sketch.tsv
      -k  12  flanking kmer length. For k=12, the total window size is 25.
      --mincount 2  Minimum number of left/right flanking kmer combinations
      
Second, find distances

      Usage: comboverDist.pl kmers1.tsv kmers2.tsv
