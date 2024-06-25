#!/usr/bin/perl -w

use strict;

#Usage: $0 gene_file

my %genesnp;
open (IN, "$ARGV[0]");
while (<IN>)
{
    chomp;
    my @line=split(/\s+/);
    $line[0]=~s/\s+/_/g;
    push (@{$genesnp{$line[0]}}, $line[1]);
}
close (IN);

foreach my $gene ( keys %genesnp )
{
    open (OUT, ">$gene.gene");
    for (my $i=0; $i<=$#{$genesnp{$gene}}; $i++)
    {
        print OUT "$genesnp{$gene}->[$i]\n";
    }
    close (OUT);
    wait;
}
