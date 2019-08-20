#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper;

exit(main());

sub main{
  my $settings = {};
  GetOptions($settings, qw(help k mincount|min-count=i)) or die $!;
  $$settings{k} ||= 12;
  $$settings{mincount} ||= 2;

  die usage() if(!@ARGV || $$settings{help});
  
  my $infile = $ARGV[0];

  my $snp = countKmers($infile, $settings);

  # Sorted output by left/right flank
  my @sortedFlank1 = sort {$a cmp $b} keys(%$snp);
  for my $flank1(@sortedFlank1){
    my @sortedFlank2 = sort {$a cmp $b} keys(%{ $$snp{$flank1} });
    for my $flank2(@sortedFlank2){
      print join("\t", $flank1, $flank2, $$snp{$flank1}{$flank2})."\n";
    }
  }
  return 0;
  
}

sub countKmers{
  my($seqfile, $settings)=@_;

  my %snp;

  my $flankingKmerLength = $$settings{k};
  my $minCount = $$settings{mincount};

  my $windowSize = $$settings{k} * 2 + 1;

  open(my $fastqFh, "zcat \Q$seqfile\E | ") or die "ERROR zcatting $seqfile: $!";
  while(my $id = <$fastqFh>){ # burn the read ID line
    my $seq=<$fastqFh>;
    chomp($seq);
    my $length = length($seq);

    my $startPos = 0;
    my $flank2   = "";
    do{
      my $flank1   = substr($seq, $startPos, $flankingKmerLength);
      my $middleNt = substr($seq, $startPos + $flankingKmerLength, 1);
         $flank2   = substr($seq, $startPos + $flankingKmerLength + 1, $flankingKmerLength);
      $startPos += $flankingKmerLength + $flankingKmerLength + 1;

      if(length($flank2) == $flankingKmerLength){
        $snp{$flank1}{$flank2}{$middleNt}++;
      }
    } while($startPos + $windowSize < $length);

    # Burn the quality score lines
    <$fastqFh>;
    <$fastqFh>;
  }
  close $fastqFh;

  # Filter down the SNPs, removing low counts or minority counts
  my %filteredSnp;
  while(my($flank1, $flank2s) = each(%snp)){
    while(my($flank2, $ntCounts) = each(%$flank2s)){

      # For each middle nt, remove low counts and
      # get majority rules.
      while(my($nt, $count) = each(%$ntCounts)){
        if($count < $minCount){
          delete($snp{$flank1}{$flank2}{$nt});
        }
      }
      # Majority rules
      if(keys(%{ $snp{$flank1}{$flank2} }) > 0){
        my $maxCount = 0;
        my $maxNt    = "";
        while(my($nt, $count) = each(%{ $snp{$flank1}{$flank2} })){
          $filteredSnp{$flank1}{$flank2} = $nt;
          last; # TODO do actual majority rules
        }
      }
    }
  }

  return \%filteredSnp;
}


sub usage{
  "Combover: sketches a fastq file for split kmers and variant
  $0 -k 12 file.fastq.gz > sketch.tsv
  -k  12  flanking kmer length. For k=12, the total window size is 25.
  --mincount 2  Minimum number of left/right flanking kmer combinations
                before recording the variant.
  "
}
