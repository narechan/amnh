#!/usr/bin/perl -w

# -f is the fasta alignment file
# -o is the outdir
# -t is the treefile
# -c is your codeml ctl file minus the IO params
#    (seqfile, treefile, and outfile)

## TODO: still need to make this parallel 

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SeqIO;

my ($fasta, $outdir, $treefile, $configfile);
GetOptions(
	   'f|fasta=s'  => \$fasta,
	   'o|outdir=s' => \$outdir,
	   't|treefile=s' => \$treefile,
	   'c|config=s' => \$configfile,
	   );

# define stop codons (hard coded so change here)
my $stops = {};
$stops->{'TAG'} = 1;
$stops->{'TAA'} = 1;
$stops->{'TGA'} = 1;

#####MAIN######

`mkdir -p $outdir/ctlfiles`;
`mkdir -p $outdir/results`;

# get the fasta name
my $fasta_name;
if ($fasta =~m/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}
else {
    $fasta_name = $fasta;
}

print STDERR "Working on $fasta_name\n";

# get the alignment
my $alnin = Bio::AlignIO->new(-file   => "$fasta",
			      -format => 'fasta');
my $alnobj = $alnin->next_aln();

# create new alignment without stop codons
open (A, ">$outdir/ctlfiles/$fasta_name");
foreach my $seq ($alnobj->each_seq){
    my $id        = $seq->display_id;
    my $sequence  = $seq->seq;
    
    # blow apart the sequence into codons
    # and search and replace stop codons, if any
    my $newseq;
    my @codons   = $sequence =~m/(...)/g;
    foreach my $codon (@codons){
	$codon =~tr/a-z/A-Z/;
	if (exists ($stops->{$codon})){
	    $newseq .= "???";
	}
	else {
	    $newseq .= $codon;
	}
    }
    
    # write new alignment file with stop codons removed
    print A ">$id\n$newseq\n";
}
close (A);

# convert the alignment to phylip
my $alninnew = Bio::AlignIO->new(-file   => "$outdir/ctlfiles/$fasta_name",
				 -format => 'fasta');
my $alnout   = Bio::AlignIO->new(-file   => ">$outdir/ctlfiles/$fasta_name.phy",
				 -interleaved => 0,
				 -format => 'phylip');
while (my $aln = $alninnew->next_aln){
    $alnout->write_aln($aln);
}

# create control file for this job
my @controlfile;
open (C, "$configfile");
while (my $line = <C>){
    chomp $line;
    push (@controlfile, $line);
}
close (C);

my $seqfileline  = "seqfile = " . $outdir . "/ctlfiles/" . $fasta_name . ".phy";
my $treefileline = "treefile = " . $treefile;
my $outfileline  = "outfile = " . $fasta_name . ".out";
push (@controlfile, $seqfileline, $treefileline, $outfileline);

open (CN, ">$outdir/ctlfiles/$fasta_name.ctl");
foreach my $line (@controlfile){
    print CN "$line\n";
}
close (CN);

# run codeml
`codeml $outdir/ctlfiles/$fasta_name.ctl`;
`mv $fasta_name.out $outdir/results`;
`rm 2NG.* lnf rub rst*`;
