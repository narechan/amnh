#!/usr/bin/perl
    
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $itlog, $k, $outdir, $jarfile);
GetOptions(
    'h|help'         => \$help,
    'i|itlog=s'      => \$itlog,
    'k|k=s'     => \$k,
    'o|outdir=s'     => \$outdir,
    'j|farfile=s'    => \$jarfile,
           ) or pod2usage;

pod2usage if $help;

`mkdir -p $outdir`;

# get the log file name
my $itlog_name;
if ($itlog =~/\//g){
    $itlog =~m/.*\/(.*)$/;
    $itlog_name = $1;
}
else {
    $itlog_name = $itlog;
}



# run the elki optics program
`java -XX:+UseSerialGC -cp $jarfile de.lmu.ifi.dbs.elki.application.KDDCLIApplication -dbc.in $itlog -resulthandler ResultWriter -algorithm clustering.kmeans.KMeansLloyd -kmeans.k $k > $outdir/$itlog_name.kmeans`;

# parse the optics clusters
my $clusters = {};
my $idseen = {};
my $counter = -1;
open (F, "$outdir/$itlog_name.kmeans");
while (my $line = <F>){
    chomp $line;
    
    if ($line =~m/#\sCluster\:/){
	$counter++;
    }
    elsif ($line =~m/^ID/){
	my @line = split (/\s/, $line);
	if (exists ($idseen->{$line[3]})){
	    last;
	}
	else{
	    $clusters->{$counter}->{$line[3]} = 1;
	    $idseen->{$line[3]} = 1;
	}
    }
    else {
	next;
    }
}
close (F);

foreach my $clust (sort {$a <=> $b} keys %$clusters){
    foreach my $id (sort keys %{$clusters->{$clust}}){
	print "$clust\t$id\n";
    }
}
