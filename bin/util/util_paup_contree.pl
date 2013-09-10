#!/usr/bin/perl

use Getopt::Std;
use Cfi;

my %opts = ();
getopts ('m:t:o:n:', \%opts);
my $matrix = $opts{'m'};
my $outdir = $opts{'o'};
my $treedir   = $opts{'t'};
my $name   = $opts{'n'};

`mkdir -p $outdir/consensus`;

# get all the trees                                                                                         
opendir (TREES, "$outdir/trees");
my @trees = sort (readdir (TREES));
shift @trees;
shift @trees;
closedir (TREES);

my $cfiobj = Cfi->new;
$cfiobj->generate_consensus ("$outdir", $matrix, \@trees, $name);
`paup -n $outdir/consensus/all.nex`;

#open (N, ">$outdir/$matrixname.cmds");
#print N "#nexus\n";
#print N "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
#print N "execute $matrix;\n";
#print N "gettrees mode=3 file=$tree;\n";
#print N "contree all/strict=yes treefile=$outdir/$matrixname.contree;\n";
#print N "quit /warnTsave=no;\n";
#close (N);

#`paup -n $outdir/$matrixname.cmds`;

#contree all/strict=no majrule=yes usetreewts=yes treefile=bootMajRule.tre;
