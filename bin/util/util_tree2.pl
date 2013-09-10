#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('t:', \%opts);
my $treedir  = $opts{'t'};

opendir (D, "$treedir");
my @subdirs = sort (readdir (D));
shift @subdirs;
shift @subdirs;
closedir (D);

my $uniqtrees = {};
foreach my $dir (@subdirs){
    my ($q, $r) = split (/\-/, $dir);
    
    opendir (E, "$treedir/$dir");
    my @trees = sort (readdir (E));
    shift @trees;
    shift @trees;
    closedir (E);
    
    foreach my $tree (@trees){
	open (F, "$treedir/$dir/$tree");
	while (my $line = <F>){
	    chomp $line;
	    if ($line =~m/PAUP_1/){
		$uniqtrees->{$q}->{$r}->{$line}++;
	    }
	    else{
		next;
	    }
	}
	close (F);
    }
}

foreach my $query (sort {$a <=> $b} keys %$uniqtrees){
    my @string;
    foreach my $ref (sort {$a <=> $b} keys %{$uniqtrees->{$query}}){
	my $uniqs = keys %{$uniqtrees->{$query}->{$ref}};
	push (@string, $uniqs);
    }
    my $string = join "\t", @string;
    print "$query\t$string\n";
}
