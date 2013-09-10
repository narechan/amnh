#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('a:b:o:', \%opts);
my $afile  = $opts{'a'};
my $bdir  = $opts{'b'};
my $outdir = $opts{'o'};
`mkdir -p $outdir`;

my $trees = {};
opendir (D, "$bdir");
my @trees = sort (readdir(D));
shift @trees;
shift @trees;
foreach my $tree (@trees){
    my @name = split (/\./, $tree);
    
    open (T, "$bdir/$tree");
    while (my $line = <T>){
	chomp $line;
	$trees->{$name[1]} = $line;
    }
    close (T);
}
closedir (D);

$/ = "\n%%";
open (A, "$afile");
while (my $chunk = <A>){
    $chunk =~s/\%//g;
    my @lines = split (/\n/, $chunk);
    my $id = shift @lines;
    $id =~s/\-//g;
    my $ntax = @lines;
    my $matrix = join "\n", @lines;
    
    my $filename = "OID_" . $id . "_chrm.nex";
    open (O, ">$outdir/$filename");
    print O "#NEXUS\n";
    print O "BEGIN DATA;\n";
    print O "dimensions ntax=$ntax nchar=1;\n";
    print O "format missing=?\n";
    print O "symbols=\"AX\";\n";
    print O "matrix\n$matrix\n";
    print O ";\n";
    print O "end;\n";
    print O "begin trees;\n";
    print O "tree RaxML1 = $trees->{$id}\n";
    print O "end;\n";
    close (O);
}
close (A);
