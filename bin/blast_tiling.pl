#!/usr/bin/perl -w

=head1 NAME

blast_tiling.pl

=head1 SYNOPSIS

  blast_tiling.pl -- 
              

Options:

 --help        Show brief help and exit
 --blast       Is your blast result from pairwise blast (bl2seq)
 --outdir      Is your output dir

=head1 DESCRIPTION

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Search::Tiling::MapTiling;

my ($help, $blastout, $outdir);
GetOptions(
    'h|help'          => \$help,
    'o|outdir=s'      => \$outdir,
    'b|blastout=s'    => \$blastout,
    ) or pod2usage;

pod2usage if $help;

for my $option ($outdir, $blastout){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

my $in = Bio::SearchIO->new(-file   => "$blastout",
			    -format => "blast");

# only one result and one hit
my $result = $in->next_result;
my $hit    = $result->next_hit;

# get the tiling
my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
#my @qalns = $tiling->get_tiled_alns('query');
#my @salns = $tiling->get_tiled_alns('hit');

my $aln = ($tiling->get_tiled_alns)[0];
my $qseq = $aln->get_seq_by_id('query');
my $hseq = $aln->get_seq_by_id('subject');

my $qsequence = $qseq->seq;
my $hsequence = $hseq->seq;
print ">query\n$qsequence\n";
print ">subject\n$hsequence\n";



=head
# here's the concatenation:
my $qconcat_seq_obj = $qalns[0]->get_seq_by_id('query');
#my $sconcat_seq_obj = $salns[0]->get_seq_by_id('subject');

# get the originals
my $aln = ($tiling->get_tiled_alns)[0];
my $qseq = $aln->get_seq_by_id('query');
my $hseq = $aln->get_seq_by_id('subject');
foreach my $feat ($qseq->get_SeqFeatures) {
    my $org_start = ($feat->get_tag_values('query_start'))[0];
    my $org_end = ($feat->get_tag_values('query_end'))[0];
    
    # original fragment as represented in the tiled alignment:
    my $org_fragment = $feat->seq->seq;
    print "$org_start\t$org_end\t$org_fragment\n";
}

foreach my $feat ($hseq->get_SeqFeatures) {
    my $org_start = ($feat->get_tag_values('subject_start'))[0];
    my $org_end = ($feat->get_tag_values('subject_end'))[0];
    
    # original fragment as represented in the tiled alignment:
    my $org_fragment = $feat->seq->seq;
    print "$org_start\t$org_end\t$org_fragment\n"
}
=cut
