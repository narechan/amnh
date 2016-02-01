#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('d:r:', \%opts);
my $dir  = $opts{'d'};
my $ref  = $opts{'r'};

opendir (D, "$dir");
my @res = sort(readdir (D));
shift @res;
shift @res;
closedir (D);

my $strings = {};
foreach my $resfile (@res){
    my $data = {};
    open (F, "$dir/$resfile");
    while (my $line = <F>){
	chomp $line;
	my @line = split (/\,/, $line);
	$data->{$line[1]} = $line[6];
    }
    my @genes;
    my @values;
    my @ref;
    foreach my $key (sort keys %$data){
	push (@genes, $key);
	push (@values, $data->{$key});
	push (@ref, 0);
    }
    my $genes = join "\t", @genes;
    my $values = join "\t", @values;
    my $refs = join "\t", @ref;
    $strings->{'genes'} = $genes;
    $strings->{$ref} = $refs;
    $strings->{$resfile} = $values;
}

foreach my $string (keys %$strings){
    print "$string\t$strings->{$string}\n";
}

