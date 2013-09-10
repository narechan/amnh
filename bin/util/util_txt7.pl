#!/usr/bin/perl

# field4 = P1
# field6 = P2
# field8 = P3
# field11 = Pholm

use Getopt::Std;

my %opts = ();
getopts ('p:i:f:t:l:', \%opts);
my $pval  = $opts{'p'};
my $indir    = $opts{'i'};
my $field = $opts{'f'};
my $tuber = $opts{'t'};
my $lookup = $opts{'l'};

my $class = {};
my $func = {};
open (T, "$tuber");
while (my $line = <T>){
    chomp $line;
    my @line = split (/\t/, $line);
    $class->{$line[1]} = $line[0];
    $func->{$line[1]} = $line[3];
}
close (T);

my $lin = {};
open (L, "$lookup");
while (my $line = <L>){
    chomp $line;
    my @line = split (/\t/, $line, 2);
    $lin->{$line[0]} = $line[1];
}
close (L);


opendir (D, "$indir");
my @exptfiles = sort (readdir(D));
shift @exptfiles;
shift @exptfiles;
closedir (E);

my $fieldarray = $field - 1;

foreach my $file (@exptfiles){
    my @file = split (/\./, $file);
    my $filename = shift (@file);
    
    open (F, "$indir/$file");
    while (my $line = <F>){
	chomp $line;
	next if $line =~m/^Branch/;
	
	my @data = split (/\,/, $line);
	my $data = join "\t", @data;
	(print "$filename\t$lin->{$data[0]}\t$func->{$filename}\t$class->{$filename}\t$data\n") if ($data[$fieldarray] <= $pval);
    }
    close (F);
}
