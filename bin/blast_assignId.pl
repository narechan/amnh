#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::AlignIO;
use Bio::Index::Fasta;

my ($help, $queryfasta, $reffasta, $config);
GetOptions(
    'q|query=s'      => \$queryfasta,
    'r|ref=s'        => \$reffasta,
    'c|config=s'     => \$config,
    ) or pod2usage;

pod2usage if $help;

# NOTE that a blastdb must exist for the reffasta in the same
# path as the fasta file

# if using tblastn, the ref should be the nucs and the query
# should be the aas

#####MAIN#####

# get the config                                                                                      
my $conf = parse_config ($config);

# create hash lookup of the reffasta
# needed for when blast does not index the query deflines 
my $reffastalookup = {};
my $reffasta_name;
if ($reffasta =~/\//g){
    $reffasta =~m/.*\/(.*)$/;
    $reffasta_name = $1;
}
else {
    $reffasta_name = $reffasta;
}

my $rcounter;
my $seqin0 = Bio::SeqIO->new (-format=>'Fasta', -file=>"$reffasta");
while (my $sequence_obj = $seqin0->next_seq()){
    $rcounter++;
    my $seq      = $sequence_obj->seq();
    my $id = $rcounter . "_" . $reffasta_name;
    $reffastalookup->{$id} = $seq;
}

# create sequence database for reffasta
my $index = Bio::Index::Fasta->new(-filename => $reffasta . ".idx", -write_flag => 1);
$index->make_index($reffasta);

# do the blast for each query
my $newfasta = {};
my $seqlookup = {};
my $qcounter = 0;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$queryfasta");
while (my $sequence_obj = $seqin->next_seq()){
    $qcounter++;
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
#    print STDERR "$qcounter\tWorking on $id\n";
    
    # print out a temp fasta file
    `mkdir -p /tmp/blast/$qcounter`;
    open (O, ">/tmp/blast/query.$qcounter");
    print O ">$id\n$seq\n";
    close (O);

    # blast
    my $name = blast ($conf, "/tmp/blast/query.$qcounter", "$reffasta", "/tmp/blast/$qcounter");

    # parse
    my $hitkeep;
    my $fracid;
    my $in = Bio::SearchIO->new(-file   => "/tmp/blast/$qcounter/$name",
				-format => "blast");
    while (my $result = $in->next_result){
	my $qname = $result->query_name;
	
	if ($result->hits){ 
	    
	    # assume the top hit is the ID correspondence we want
	    while (my $hit = $result->next_hit){
		my $hname = $hit->name;
		my $evalue = $hit->significance;
		
		$hitkeep = $hname;
		
		# get the fraction of identities in the first HSP
		# (and only HSP if clean, perfect hit)
		while (my $hsp = $hit->next_hsp){
		    $fracid = $hsp->frac_identical;
		    last;
		}
		last;
	    }

	
	    print STDERR "$qcounter\tWorking on $id\t$fracid\n";
	    
	    # get the refseq and storex
	    my $hseq = $index->fetch($hitkeep);
	    if ($hseq){
		$newfasta->{$id} = $hseq->seq;
	    }
	    else {
		$newfasta->{$id} = $reffastalookup->{$hitkeep};
	    }
	}
	else {
	    print STDERR "$qcounter\tWorking on $id\tNO HITS\n";
	    $newfasta->{$id} = "";
	}
    }
}
    
# print out the fasta file with reassigned ids
foreach my $nid (keys %$newfasta){
    print ">$nid\n$newfasta->{$nid}\n";
}

`rm -rf /tmp/blast`;
    
#####SUBS#####
sub blast{
    my $conf   = shift;
    my $fasta  = shift;
    my $db     = shift;
    my $out    = shift;

    my $fasta_name;
    if ($fasta =~/\//g){
        $fasta =~m/.*\/(.*)$/;
        $fasta_name = $1;
    }

    else {
        $fasta_name = $fasta;
    }

    my $db_name;
    if ($db =~/\//g){
        $db =~m/.*\/(.*)$/;
        $db_name = $1;
    }

    else {
        $db_name = $db;
    }

    my $cmd = $conf->{BLASTALL};
    $cmd .= " $conf->{PARAMETERS}";
    $cmd .= " -i $fasta";
    $cmd .= " -d $db";
    $cmd .= " -o $out/$fasta_name.$db_name.out";
    `$cmd`;

    return ("$fasta_name.$db_name.out");
}

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }

    # check for blastall in user's path if not in config                                                   
    if (! exists($config{"BLASTALL"})){
        my $blastall = `which blastall`;
        chomp $blastall;

        if (defined ($blastall)){
            $config{"BLASTALL"} = $blastall;
        }
        else {
            die "no BLASTALL specified!\n";
        }
    }

    return (\%config);
    close (F);
}
