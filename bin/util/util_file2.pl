#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:d:', \%opts);
my $dir  = $opts{'d'};
my $list = $opts{'l'};

opendir (D, "$dir");
my @files = sort(readdir (D));
shift @files;
shift @files;
closedir (D);


my $hash = {};
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    my ($acc, $name) = split (/\t/, $line);
#    $hash->{$acc} = $name;
    $hash->{$name} = $acc;
}
close (L);


foreach my $file (@files){
    `mv $dir/$file $dir/$hash->{$file}`;
#    open (F, "$dir");
#    while (my $line = <F>){
#	chomp $line;
#	my ($acc, $seq) = split (/\t/, $line);
#	print "$acc\n$seq\n";
#    }
#    close (F);
}

