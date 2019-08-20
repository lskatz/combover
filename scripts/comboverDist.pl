#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help} || @ARGV < 2);

  distance(\@ARGV, $settings);
  
  return 0;
}

sub distance{
  my($files, $settings) = @_;
  
  my @dist;
  for my $filename(@$files){
    my $dist = readDistFile($filename, $settings);
    push(@dist, $dist);
  }

  for(my $i=0;$i<@dist;$i++){
    for(my $j=0;$j<@dist;$j++){
      my %jaccard = compareDistances($dist[$i], $dist[$j], $settings);
      $jaccard{jaccard} = sprintf("%0.2e", $jaccard{jaccard});
      print join("\t", $$files[$i], $$files[$j], $jaccard{jaccard}, $jaccard{diff}, $jaccard{total})."\n";
    }
  }
}

sub readDistFile{
  my($filename, $settings) = @_;
  
  my %dist;

  open(my $fh, $filename) or die "ERROR: could not open $filename: $!";
  while(<$fh>){
    chomp;
    my($kmer1,$kmer2,$nt) = split(/\t/,$_);
    $dist{$kmer1}{$kmer2}=$nt;
  }
  close $fh;
  return \%dist;
}

sub compareDistances{
  my($dist1,$dist2,$settings)=@_;

  my $diff = 0;
  my $total= 0;
  # Compare kmer set #1 to set #2
  while(my($kmer1,$kmer2s)=each(%$dist1)){
    while(my($kmer2,$nt) = each(%$kmer2s)){
      if(defined($$dist2{$kmer1}{$kmer2})){
        $total++;
        if($$dist1{$kmer1}{$kmer2} ne $$dist2{$kmer1}{$kmer2}){
          $diff++;
        }
      }
    }
  }
  # Compare kmer set #2 to set #1
  while(my($kmer1,$kmer2s)=each(%$dist2)){
    while(my($kmer2,$nt) = each(%$kmer2s)){
      if(defined($$dist1{$kmer1}{$kmer2})){
        $total++;
        if($$dist1{$kmer1}{$kmer2} ne $$dist2{$kmer1}{$kmer2}){
          $diff++;
        }
      }
    }
  }
  
  my $jaccard = $diff/$total;
  return (jaccard=>$jaccard, diff=>$diff, total=>$total) if wantarray;
  return $jaccard;
}

sub usage{
  "$0: distances between two split kmer sets
  Usage: $0 kmers1.tsv kmers2.tsv
  --help   This useful help menu
  "
}
