#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;

$ENV{PATH}="./scripts:$ENV{PATH}";

for my $file(glob("t/data/*.fastq.gz")){
  system("comboverSketch.pl $file > $file.tsv");
  die "ERROR sketching $file" if $?;
}

my @dist = `comboverDist.pl t/data/*.fastq.gz.tsv`;
die "ERROR with distances" if $?;
chomp(@dist);

subtest "sample1/sample1" => sub{
  my @line1 = split(/\t/, $dist[0]);
  is($line1[3], 0, "zero distance self/self");
  is($line1[4], 20498, "20498 kmers self/self");
};

subtest "sample1/sample2" => sub{
  my @line2 = split(/\t/, $dist[1]);
  is($line2[3], 8, "8 differences between sample1/sample2");
  is($line2[4], 1754, "1754 kmers b/n sample1/sample2");
};
