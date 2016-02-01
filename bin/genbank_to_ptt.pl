#!/usr/bin/perl -w

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;

my ($help, $gb, $type);
GetOptions(
    'h|help'          => \$help,
    'b|gb=s'          => \$gb,
    't|type=s'        => \$type, # CDS or RNA
	   ) or pod2usage;
pod2usage if $help;

for my $option ($gb, $type){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

my $seqin = Bio::SeqIO->new(-file=>$gb, -format=>'Genbank');       
my $seqobj = $seqin->next_seq;
my $chrom  = $seqobj->display_id;

print $chrom, " - 1..", $seqobj->length,"\n";
my @cds;
if ($type eq "CDS"){
    @cds = grep { $_->primary_tag eq $type } $seqobj->get_SeqFeatures;
    print scalar(@cds)," proteins\n";
}
elsif ($type eq "RNA"){
    my @rrna = grep { $_->primary_tag eq 'rRNA' } $seqobj->get_SeqFeatures;
    my @trna = grep { $_->primary_tag eq 'tRNA' } $seqobj->get_SeqFeatures;
    push (@cds, @rrna);
    push (@cds, @trna);
    print scalar(@cds)," RNAs\n";
}
else {
    print STDERR "Unrecognized type\n";
    die;
}

print join("\t", qw(Location Strand Length PID Gene Synonym Code COG Product)),"\n";

# cycle through the features
my $counter = 0;
foreach my $feat (sort {$a->start <=> $b->start} @cds){
    $counter++;
#    print STDERR "$counter\n";
    
    # parse the data                                                          
    my $name;
    if ($feat->has_tag("db_xref")){
        my @names = $feat->get_tag_values("db_xref");
	foreach my $n (@names){
	    ($name = $n) if ($n =~m/^SEED/);
	}
    }
    else {
        $name = "NONAME_$counter";
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
    my $sestring = $start . ".." . $end;
    my $strand = $feat->strand;
    my $len = ($feat->length / 3) - 1;
    print "$sestring\t$strand\t$len\t-\t$name\t$name\t-\t-\t$desc\n";
}
