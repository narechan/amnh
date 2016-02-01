#!/usr/bin/perl

use Getopt::Std;
use Parallel::ForkManager;

my %opts = ();
getopts ('m:o:', \%opts);
my $matrix = $opts{'m'};
my $outdir = $opts{'o'};

my $matrix_name;
if ($matrix =~/\//g){
    $matrix =~m/.*\/(.*)$/;
    $matrix_name = $1;
}

else {
    $matrix_name = $fasta;
}


open (N, ">$outdir/$matrix_name.cmds");
print N "#nexus\n";
print N "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
print N "set maxtrees=200 increase=no;\n";
print N "log start file=$outdir/$matrix_name.out replace=yes;\n";
print N "execute $matrix;\n";
print N "exclude all;\n";
print N "include all;\n";
print N "hsearch swap=tbr addseq=random nreps=100 multrees=no;\n";
print N "savetree replace=yes file=$outdir/$matrix_name.tre;\n";
print N "log stop;\n";
print N "quit /warnTsave=no;\n";
close (N);

`paup -n $outdir/$matrix_name.cmds`;
