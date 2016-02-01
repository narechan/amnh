#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('d:', \%opts);
my $dir  = $opts{'d'};

opendir (D, "$dir");
my @files = sort(readdir (D));
shift @files;
shift @files;
closedir (D);

my $data = {};
my $queries = {};
foreach my $file (@files){
    open (F, "$dir/$file");
    while (my $line = <F>){
	chomp $line;
	my @line = split (/\t/, $line);
	if ($line[5] == 1){
	    $data->{$file}->{$line[0]} = $line[4];
	    $queries->{$line[0]}++;
	}
	else {
	    next;
	}
    }
    close (F);
}

my @header;
foreach my $query (sort keys %$queries){
    push (@header, $query);
}
my $header = join "\t", @header;
print "\t$header\n";

foreach my $strain (sort keys %$data){
    my @string;
    foreach my $query (sort keys %$queries){
	push (@string, $data->{$strain}->{$query});
    }
    my $string = join "\t", @string;
    print "$strain\t$string\n";
}


###subs###

sub file_count{
    my $file = shift;

    my $counter = 0;
    open (F, "$file");
    while (my $line = <F>){
	chomp $line;
	$counter++;
    }
    close (F);
    
    return ($counter);
}
