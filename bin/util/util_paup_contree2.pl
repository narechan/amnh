#!/usr/bin/perl

use Getopt::Std;
use Cfi;

my %opts = ();
getopts ('m:t:o:', \%opts);
my $matrix = $opts{'m'};
my $outdir = $opts{'o'};
my $tetree   = $opts{'t'};

`mkdir -p $outdir`;

open (N, ">tmp.cmds");
print N "#nexus\n";
print N "log start file=log replace=yes;\n";
print N "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
print N "execute $matrix;\n";
print N "gettrees mode=7 file=$tetree;\n";
print N "contree all/ indices=yes;\n";
print N "quit /warnTsave=no;\n";
close (N);

`paup -n tmp.cmds`;
#`rm tmp.cmds`;

#contree all/strict=no majrule=yes usetreewts=yes treefile=bootMajRule.tre;
