#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:o:', \%opts);
my $indir = $opts{'i'};
my $outdir = $opts{'o'};

opendir (D, "$indir");
my @files = sort (readdir (D));
shift @files;
shift @files;
closedir (D);

`mkdir -p $outdir`;

foreach my $file (@files){
    
    open (O, ">$outdir/$file.tre");
    
    open (F, "$indir/$file");
    my $tree;
    while (my $line = <F>){
	chomp $line;
	$tree = $line;
    }
    close (F);
    
    print O "#NEXUS\n";
    print O "Begin trees;\n";
    print O "tree $file=$tree\n";
    print O "End;\n";

    close (O);
}
