#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('d:o:', \%opts);
my $dir  = $opts{'d'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

opendir (D, "$dir");
my @resdirs = sort(readdir (D));
shift @resdirs;
shift @resdirs;
closedir (D);

my %genomes_encountered;
foreach my $resdir (@resdirs){
    my @data = split (/\-/, $resdir);
    my $genome = $data[0];
    my $tech = $data[3];

    # technology filter                                                                               
    next unless ($tech eq "solexa");

    # genome filter
    if (exists ($genomes_encountered{$genome})){
	next;
    }
    else {
	$genomes_encountered{$genome} = 1;
    }

    print STDERR "$genome\n";
#    `mkdir $outdir/$genome`;
    
    `tar xzf $dir/$resdir/metasim.tar.gz`;
    `mv results $dir/$resdir`;

    opendir (R, "$dir/$resdir/results/$resdir/metasim/reads");
    my @readfiles = sort {$a <=> $b} (readdir (R));
    shift @readfiles;
    shift @readfiles;
    closedir (R);

    `cp $dir/$resdir/results/$resdir/metasim/reads/$readfiles[5] $outdir`;
    `mv $outdir/$readfiles[5] $outdir/$genome.1X`;
    `cp $dir/$resdir/results/$resdir/metasim/reads/$readfiles[13] $outdir`;
    `mv $outdir/$readfiles[13] $outdir/$genome.5X`;
    `cp $dir/$resdir/results/$resdir/metasim/reads/$readfiles[19] $outdir`;
    `mv $outdir/$readfiles[19] $outdir/$genome.10X`;
    `cp $dir/$resdir/results/$resdir/metasim/reads/$readfiles[21] $outdir`;
    `mv $outdir/$readfiles[21] $outdir/$genome.20X`;

    `rm -rf $dir/$resdir/results`;
    
}
