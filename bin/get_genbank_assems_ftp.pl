#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;
use Bio::SeqIO;

my ($help, $list, $outdir); #, $url);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   'o|outdir=s'   => \$outdir,
           'u|url=s'      => \$url,
	   );

pod2usage(2) if $help;

`mkdir -p $outdir/reports`;
`mkdir -p $outdir/assems`;
`mkdir -p $outdir/proteins`;
`mkdir -p $outdir/gffs`;

my $key = {};
my $names = {};
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    print STDERR "Working on $line\n";
    
    # get the assembly report
    my $rprtcmd = "curl $url/$line/" . "$line" . "_assembly_report.txt >  $outdir/reports/$line.report";
#    my $rprtcmd = "curl ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/$line.assembly.txt > $outdir/reports/$line.report";
    `$rprtcmd`;
    
    # parse the assmbly report and get the org name
    my $name;
    my $assemblyname;
    my $stuff;
    my $assemstuff;
    open (A, "$outdir/reports/$line.report");
    while (my $line = <A>){
        chomp $line;
        if ($line =~m/Organism name\:/){
	    ($stuff, $name) = split (/\:/, $line);
	    $name =~s/[[:punct:]]/_/g;
	    $name =~s/^\s+//g;
	    $name =~s/\s+$//g;
	    $name =~s/\s/_/g;
	}
	elsif ($line =~m/Assembly Name\:/){
	    ($assemstuff, $assemblyname) = split (/\:/, $line);
            $assemblyname =~s/[[:punct:]]/_/g;
            $assemblyname =~s/^\s+//g;
            $assemblyname =~s/\s+$//g;
            $assemblyname =~s/\s/_/g;
        }

	else {
	    next;
	}
    }
    close (A);
    
    if (exists ($names->{$name})){
	$names->{$name}++;
	$name .= $names->{$name};
    }
    else{
	$names->{$name} = 1;
    }

    $key->{$line} = $name;
    
    # get the contigs
    my $ctgcmd = "curl $url/$line/" . "$line" . "_genomic.fna.gz | gzip -d > $outdir/assems/$name.fa";
#    my $string = "$line" . "_" . "$assemblyname";
#    my $ctgcmd = "curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$string/$string" . "_genomic.fna.gz | gzip -d > $outdir/assems/$name.fa";
    `$ctgcmd`;
    
    # get the gffs
    my $gffcmd = "curl $url/$line/" . "$line" . "_genomic.gff.gz | gzip -d > $outdir/gffs/$name.gff";
    `$gffcmd`;
    
    # get the proteins
    my $ptcmd = "curl $url/$line/" . "$line" . "_protein.faa.gz | gzip -d > $outdir/proteins/$name.fa";
    `$ptcmd`;
}
close (L);

open (K, ">$outdir/key");
foreach my $a (sort %$key){
    print K "$a\t$key->{$a}\n";
}
close (K);
