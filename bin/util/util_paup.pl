#!/usr/bin/perl

use Getopt::Std;
use Parallel::ForkManager;

my %opts = ();
getopts ('p:i:o:', \%opts);
my $indir = $opts{'i'};
my $outdir = $opts{'o'};
my $procs  = $opts{'p'};

`mkdir -p $outdir/trees`;
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;

opendir (E, "$indir");
my @files = sort (readdir(E));
shift @files;
shift @files;
closedir (E);
my $pm = Parallel::ForkManager->new($procs);
foreach my $file (@files){
    $pm->start and next;
    warn "$file\n";
    open (N, ">$outdir/cmds/$file.cmds");
    print N "#nexus\n";
    print N "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
    print N "set maxtrees=200 increase=no;\n";
    print N "log start file=$outdir/logs/$file.out replace=yes;\n";
    print N "execute $indir/$file;\n";
    print N "exclude all;\n";
    print N "include all;\n";
    print N "hsearch swap=tbr addseq=random nreps=20 multrees=no;\n";
    print N "savetree replace=yes file=$outdir/trees/$file.tre;\n";
    print N "log stop;\n";
    print N "quit /warnTsave=no;\n";
    close (N);

    `paup -n $outdir/cmds/$file.cmds`;
    $pm->finish;
}
$pm->wait_all_children;
