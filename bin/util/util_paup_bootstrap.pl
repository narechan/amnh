#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('m:o:n:', \%opts);
#getopts ('m:o:r:', \%opts);
my $matrix = $opts{'m'};
my $outdir = $opts{'o'};
my $reps   = $opts{'n'};
#my $root   = $opts{'r'};
`mkdir -p $outdir`;

my @matrix = split (/\//, $matrix);
my $matrixname = pop @matrix;

open (N, ">$outdir/$matrixname.cmds");
print N "#nexus\n";
print N "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
print N "log start file=$outdir/$matrixname.out replace=yes;\n";
print N "execute $matrix;\n";
#print N "delete TW20 ATCC_BAA_39 JKD6009 JKD6008 T0131;\n";
#print N "outgroup $root;\n";
#print N "bootstrap nreps=$reps treefile=$matrixname.tre/swap=tbr addseq=random nreps=100 multrees=no;\n";
#print N "bootstrap nreps=$reps treefile=$matrixname.tre/swap=tbr addseq=random nreps=20 multrees=no;\n";
print N "set criterion=distance;\n";
print N "bootstrap nreps=$reps treefile=$matrixname.tre search=nj/brlens=yes;\n";
print N "savetrees from=1 to=1 savebootp=nodelabels brlens=yes file=$matrixname.contre;\n";
print N "log stop;\n";
print N "quit /warnTsave=no;\n";
close (N);

`paup -n $outdir/$matrixname.cmds`;
`mv $matrixname.tre $outdir/`;
`mv $matrixname.contre $outdir/`;
