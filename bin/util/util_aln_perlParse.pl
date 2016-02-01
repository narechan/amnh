#!/usr/bin/perl -w

# -i is the alignment infile

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile);
GetOptions(
	   'i|infile=s'      => \$infile,
	   );

#####MAIN#####

my $aln = {};
open (A, "$infile");
while (my $line = <A>){
    chomp $line;
#    print STDERR "$line\n";
    next unless ($line);
    next if ($line =~m/^\[/);
    my ($acc, $snp) = split (/\t/, $line);
    push (@{$aln->{$acc}}, $snp);
}
close (A);

foreach my $tax (sort keys %$aln){
    my $seqstring = join "", @{$aln->{$tax}};
    my $len = length $seqstring;
    print "$tax\t$seqstring\n";
    print STDERR "$tax\t$len\n";
}
