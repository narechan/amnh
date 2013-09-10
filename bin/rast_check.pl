#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use Bio::SeqIO;

my ($help, $aafile, $nucfile, $gff3dir, $fastadir, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fastadir=s'       => \$fastadir,
    'o|outdir=s'      => \$outdir,
    'a|aafile=s'      => \$aafile,	   
    'n|nucfile=s'     => \$nucfile,
    'g|gff3dir=s'     => \$gff3dir,	   
    ) or pod2usage;

pod2usage if $help;

for my $option ($fastadir, $outdir, $aafile, $nucfile, $gff3dir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/pts`;
`mkdir -p $outdir/nucs`;

#####MAIN#####

# create the rast aa and nuc inidices
my $nucindex = Bio::Index::Fasta->new(-filename => $nucfile . ".idx", -write_flag => 1);
$nucindex->make_index($nucfile);
my $aaindex = Bio::Index::Fasta->new(-filename => $aafile . ".idx", -write_flag => 1);
$aaindex->make_index($aafile);

# associate contigs with accessions through gff3 file
my $gff3data = {};
opendir (G, "$gff3dir");
my @gff3dirs = sort (readdir (G));
shift @gff3dirs;
shift @gff3dirs;
closedir (G);

foreach my $gff3file (@gff3dirs){
    open (GFF, "$gff3dir/$gff3file");
    while (my $line = <GFF>){
        next if ($line =~m/^\#/);
        chomp $line;

        my ($chrom,
            $source,
            $type,
            $start,
            $end,
            $score,
            $strand,
            $phase,
            $attr) = split (/\t/, $line);

        my @attrs = split (/\;/, $attr);
	foreach my $attr (@attrs){
	    my ($key, $value) = split (/\=/, $attr);
	    if ($key eq "ID"){
		($gff3data->{$chrom}->{$value} = 1) unless ($value =~m/rna/);
	    }
	    else {
		next;
	    }
	}
    }
    close (GFF);
}

# cycle through contig files and collate pts and nucs
opendir (F, "$fastadir");
my @fastadirs = sort (readdir (F));
shift @fastadirs;
shift @fastadirs;
closedir (F);

foreach my $fastafile (@fastadirs){
    open (NOUT, ">$outdir/nucs/$fastafile");
    open (POUT, ">$outdir/pts/$fastafile");
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fastadir/$fastafile");
    while (my $sequence_obj = $seqin->next_seq()){
        my $contigid = $sequence_obj->display_id();
	print STDERR "$contigid\n";

	foreach my $pegid (keys %{$gff3data->{$contigid}}){
	    print STDERR "$pegid\n";
	    my $nucobj = $nucindex->fetch($pegid);
	    my $ptobj  = $aaindex->fetch($pegid);
	    my $nucseq = $nucobj->seq;
	    $nucseq =~tr/a-z/A-Z/;
	    my $ptseq  = $ptobj->seq;
	    $ptseq =~tr/a-z/A-Z/;

	    $pegid =~s/\|//g;
	    $pegid =~s/\.//g;
	    print NOUT ">$fastafile#$pegid\n$nucseq\n";
	    print POUT ">$fastafile#$pegid\n$ptseq\n";
	}
    }
    close (NOUT);
    close (POUT);
}
	
		
