#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $infile    = $opts{'i'};

my $charset = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    $line =~tr/a-z/A-Z/;
    my ($partition, $coords) =
	split (/\=/, $line);

    $partition =~s/charset\s*//ig;
    $partition =~s/\s+//g;

    if ($coords =~m/\-/){
	$coords =~s/\;//g;
	$coords =~s/^\s//;

	$charset->{$partition} = $coords;
    }
    else {
	next;
    }
}
close (F);

open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    $line =~tr/a-z/A-Z/;
    my ($partition, $coords) =
        split (/\=/, $line);

    $partition =~s/charset\s*//ig;
    $partition =~s/\s+//g;

    if ($coords =~m/\-/){
	print "CHARSET $partition=$coords\n";
    }
    else {
	$coords =~s/\t+\;//g;
	$coords =~s/\s+\;//g;
	$coords =~s/\;//g;
	$coords =~s/^\t+//;
	$coords =~s/^\s+//;
	my @parts = split (/\s+/, $coords);
	my @cs;
	foreach my $part (@parts){
	    if (exists ($charset->{$part})){
		push (@cs, $charset->{$part});
	    }
	    else{
		print STDERR "$part undefined\n";
	    }
	}
	my $cs = join " ", @cs;
	print "CHARSET $partition=$cs;\n";
    }
}
close (F);

