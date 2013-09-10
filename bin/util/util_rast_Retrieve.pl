#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:o:', \%opts);
my $list  = $opts{'l'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

open (F, "$list");
while (my $line = <F>){
    chomp $line;
    my ($jobid, $id, $acc) = split (/\t/, $line);
    print STDERR "Working on $id\n";
    `svr_retrieve_RAST_job anarechania kringle $jobid nucleic_acid > $outdir/$acc.fa`;

#    `svr_retrieve_RAST_job anarechania kringle $jobid rast_tarball > $outdir/$acc.tar.gz`;
#    `cp $outdir/$id/Features/peg/fasta $outdir/$acc.fa`;
}
close (F);
