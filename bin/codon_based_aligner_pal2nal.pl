#!/usr/bin/perl -w

# -a is the input pt aln dir
# -i is the input pt aln format
# -n is a fasta file containing at least all corresponding nuc sequences (may contain other stuff to)
# -l is the lookup file of OG to accesssions
# -f is the codon aln out format
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::SeqIO;
use Bio::Index::Fasta;

my ($alndir, $informat, $nucfasta, $outformat, $lookupfile, $outdir);
GetOptions(
	   'a|alndir=s'      => \$alndir,
	   'i|informat=s'    => \$informat,
	   'n|nucfasta=s'     => \$nucfasta,
	   'f|outformat=s'    => \$outformat,
	   'l|lookup=s'       => \$lookupfile,
	   'o|outdir=s'       => \$outdir,
	   );

#####MAIN######

`mkdir -p $outdir`;

# index the nuc fasta
my $index = Bio::Index::Fasta->new(-filename => $nucfasta . ".idx", -write_flag => 1);
$index->make_index($nucfasta);

# get the nuc accs
my $lookup = {};
open (L, "$lookupfile");
while (my $line = <L>){
    chomp $line;
    my ($og, $accs) = split (/\t/, $line);
    my @accs = split (/\,/, $accs);
    $lookup->{$og} = [@accs];
}
close (L);

# read in all aa alns and do codon alignment
opendir (D, "$alndir");
my @ogs = sort (readdir (D));
shift @ogs;
shift @ogs;
closedir (D);
    
foreach my $og (@ogs){
    print STDERR "$og\n";
    $og =~s/\.fa//g;
    my @accs = @{$lookup->{$og}};
    
    my $nucs = {};
    my $nucsizes = {};
    
    open (T, ">zzzz.fa");
    open (Z, ">yyyy.fa");
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$alndir/$og.fa");
    while (my $sequence_obj = $seqin->next_seq()){
        my $id       = $sequence_obj->display_id();
	my $seq      = $sequence_obj->seq();
	
	my $report;
	foreach my $acc (@accs){
	    my ($sp, $ac) = split (/\#/, $acc);
	    if ($sp eq $id){
		print T ">$acc\n$seq\n";
		$report = $acc;
	    }
	    else{
		next;
	    }
	}
	
	my $nucobj = $index->fetch($report);
	my $nid       = $nucobj->id();
	my $nseq      = $nucobj->seq();
	print Z ">$nid\n$nseq\n";
    }
    close (T);
    close (Z);

    # run pal2nal
    `~apurva/packages/pal2nal.v14/pal2nal.pl zzzz.fa yyyy.fa -output fasta > $outdir/$og`;
    `rm zzzz.fa yyyy.fa`;
}











=head
    foreach my $acc (@accs){
#	if ($acc =~m/^2/){
#	    $acc = "Sa_" . $acc;
#	}
	my $nucobj = $index->fetch($acc);
	my ($strain, $id) = split (/\#/, $acc);

	# take the largest gene for the codon alignment
	my $len = $nucobj->length;
	if (! exists ($nucsizes->{$strain})){
	    $nucs->{$strain} = $nucobj;
            $nucsizes->{$strain} = $len;
	}
	else {
	    if ($len > $nucsizes->{$strain}){
		$nucs->{$strain} = $nucobj;
		$nucsizes->{$strain} = $len;
	    }
	    else{
		next;
	    }
	}
    }

    my $aa_aln = Bio::AlignIO->new(-file   => "$alndir/$og.fa",
				   -format => "$informat");
    my $aa_obj = $aa_aln->next_aln();   
    my $dna_aln = aa_to_dna_aln($aa_obj, $nucs);
    my $alnout = Bio::AlignIO->new(-file   => ">$outdir/$og.fa",
				   -format => "$outformat");
    $alnout->write_aln($dna_aln);
}
=cut
