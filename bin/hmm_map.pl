#!/usr/bin/perl

# This program will score 
# hmms (oid based, usually),
# against a reference fasta file
# and assign coordinates based on
# a minimum cutoff

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Index::Fasta;

my %opts = ();
getopts ('h:f:c:o:d:b:g:', \%opts);
my $fasta  = $opts{'f'}; # input fasta file to score against
my $gb     = $opts{'g'}; # input genbank to score against
my $oiddir = $opts{'d'}; # directory of oid families as hmms
my $cutoff = $opts{'b'}; # using bit score cutoff here, ie 50
my $outdir = $opts{'o'}; # output
my $config = $opts{'c'}; # contains any params for hmmsearch
#my $datafile = $opts{'x'}; # datafile (right now ild, but can be anything)

`mkdir $outdir`;

# parse the config
my $conf = {};
open (F, "$config");
while (my $line = <F>){
    chomp $line;
    my ($key, $value) = split (/\=/, $line, 2);
    $conf->{$key} = $value;

}
close (F);

# process the ild data file
=head
my $data = {};
open (D, "$datafile");
while (my $line = <D>){
    chomp $line;
    my @line = split (/\t/, $line);
    $data->{$line[0]} = $line[2];
}
close (D);
=cut

print STDERR "Processing the reference\n";
# index the fasta file                                                              
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

my $seqin = Bio::SeqIO->new(-file=>$gb, -format=>'Genbank');
my $seqobj = $seqin->next_seq;
my $chrom  = $seqobj->display_id;
my $seqLen = $seqobj->length;

# cycle through the features                                                         
my $coords = {};
open (P, ">$outdir/pt.fa");
foreach my $feat ($seqobj->get_SeqFeatures){
    my $type = $feat->primary_tag;
    next unless ($feat->primary_tag =~m/CDS|tRNA|rRNA|transcript/);

    # parse the data                                                                 
    my $name;
    if ($feat->has_tag("gene")){
        my @names = $feat->get_tag_values("gene");
        $name = join "_", @names;
    }
    elsif ($feat->has_tag("locus_tag")){
        my @names = $feat->get_tag_values("locus_tag");
        $name = join "_", @names;
    }
    else {
        print STDERR "No name!\n";
        die;
    }

    my $desc;
    if ($feat->has_tag("product")){
        my @names = $feat->get_tag_values("product");
        $desc = join "_", @names;
    }
    else {
        $desc = "NODESC";
    }

    my $start = $feat->start;
    my $end   = $feat->end;
    my $strand = $feat->strand;

    my $string = $start . "-" . $end;
    $coords->{$name} = $string;
    
    # harvest sequence                                                               
    my $location = Bio::Location::Simple->new(-start  => $start,
                                              -end    => $end,
                                              -strand => $strand);
    my $sequence = $index->fetch($chrom);
    my $subseq   = $sequence->subseq($location);

    my $subseq_obj = Bio::Seq->new(-seq => $subseq,
				   -alphabet => 'dna');
    my $ptsubseq_obj = $subseq_obj->translate(-frame=>0);
    my $pt = $ptsubseq_obj->seq;
    print P ">$name\n$pt\n";
}
close (P);
close (T);

# score all the oid groups against the reference fasta
opendir (D, "$oiddir");
my @ogs = sort (readdir(D));
shift @ogs;
shift @ogs;
closedir (D);

my $map = {};
foreach my $og (@ogs){
    my $hmmsearch = "hmmsearch";
    $hmmsearch   .= " $config->{'HMMSEARCH'}";
    $hmmsearch   .= " -o $outdir/$og.out --tblout $outdir/$og.tblout $oiddir/$og $outdir/pt.fa";
    `$hmmsearch`;

    my $refseq;
    my $refscore;
    open (A, "$outdir/$og.tblout");
    while (my $line = <A>){
	chomp $line;
	next if ($line =~m/\#/);
	my @line = split (/\s+/, $line);

	if ($line[5] > $score){
	    $refseq = $line[0];
	    $refscore = $line[5];
	    my $cs = $coords->{$refseq};
	    my ($s, $e) = split (/\-/, $cs);
	    for (my $i = $s; $i <= $e; $i++){
#		$map->{$i} = $refseq;
		$map->{$i} = $og;
#		$map->{$i} = $data->{$og};
	    }
	}
	else {
	    $refseq = "NONE_SIG";
	    $refscore ="NA";
	}
	last; #just want the top hit                                             
    }
    close (A);
    print STDERR "$og\t$refseq\t$refscore\n";
#    last;
}

# fill in the rest of the map with NAs and print
for (my $j = 1; $j <= $seqLen; $j++){
    if (exists ($map->{$j})){
	print "$j\t$map->{$j}\n";
	next;
    }
    else{
#	$map->{$j} = "NA";
	print "$j\tNA\n";
    }
}


