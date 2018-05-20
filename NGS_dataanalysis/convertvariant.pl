#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
  if (/^#/) {
    print;
    next;
  }
  chomp;
  my @a = split (/\t/, $_);
  my $ref = $a[3];
  my $alt = $a[4];
  for (my $i=9; $i <= $#a; $i++) {
    if ($a[$i] eq "0/0") {
      $a[$i] = $ref;
    } elsif ($a[$i] eq "1/1") {
      $a[$i] = $alt;
    } elsif ($a[$i] eq "0/1") {
      $a[$i] = "$ref$alt";
    } else {
      $a[$i] = "-";
    }
  }
  print join "\t", @a;
  print "\n";
}
