#!/usr/bin/perl -w

=head1 NAME

bl2seq_rvs.pl

=head1 SYNOPSIS

  bl2seq -- 
              
Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta file
 --ref1        Is the first reference you want to blast against (LTR sequence)   
 --ref2        Is the reference you want to blast against (entire kRV sequence)   
 --config      Is the configuration for your blast and cd-hit
 --outdir      Is your output dir
 --intlen      Is the integration site length you want to test   
    end of the virus (amount of LTR covered)

ncbi blast executables must be in your path
as must cd-hit

The config must also specify parameters for blast
except the the query and the reference seqs
Also needs parameters for CD-HIT for clustering

BLASTALL=... (optional)
PARAMETERS=...

=head1 DESCRIPTION

Given a fasta file blast all members in pairwise fashion
and parse.

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
use Bio::SeqIO;
use Bio::AlignIO;

my ($help, $fasta, $refltr, $refrv, $config, $outdir, $intlen);
GetOptions(
	   'h|help'          => \$help,
	   'f|fasta=s'       => \$fasta,
	   'o|outdir=s'      => \$outdir,
	   'c|config=s'      => \$config,
	   'a|ref1=s'        => \$refltr,
	   'b|ref2=s'        => \$refrv,
	   'i|intlen=s'      => \$intlen,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $config, $refltr, $refrv, $intlen){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# parse the fasta file
my $sequence = {};
my $seqlens  = {};
my $sin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $sin->next_seq()){
    my $id  = $sequence_obj->display_id();
    my $seq = $sequence_obj->seq();
    my $len = $sequence_obj->length();
    $sequence->{'full'}->{$id} = $seq;
    $seqlens->{$id} = $len;
}

# do the two part blast
open (RV, ">$outdir/rv.fa");
open (FLANK, ">$outdir/flank.fa");
open (NOLTRRV, ">$outdir/rv_noltr.fa");
open (NOLTRFLANK, ">$outdir/flank_noltr.fa");
open (NOLTRHITS, ">$outdir/noltrhits.fa");
open (TOOSHORT, ">$outdir/tooshort.fa");
open (LOG, ">$outdir/log");
foreach my $id1 (sort (keys %{$sequence->{'full'}})){
    print STDERR "Working on $id1\n";
    
    # generate tmp fasta file with leading 20bp primer removed
#    `mkdir -p /tmp/$id1`;
#    open (ONE, ">/tmp/$id1/$id1.rvzfa");
#    my $string = $sequence->{'full'}->{$id1};
#    $string =~s/^.{20}//; # hardcoded primer length
#    $sequence->{'noprimer'}->{$id1} = $string;
#    print ONE ">$id1\n$string\n";
#    close (ONE);
	
    # generate temp fasta file with primer (part of LTR intact)
    `mkdir -p /tmp/$id1`;
    open (ONE, ">/tmp/$id1/$id1.rvzfa");
    my $string = $sequence->{'full'}->{$id1};
    print ONE ">$id1\n$string\n";
    close (ONE);

    # ltr bl2seq
    my $ltrfilename = blast ($conf, "/tmp/$id1/$id1.rvzfa", $refltr, "$outdir/raw");
    my $ltrin = Bio::SearchIO->new(-file   => "$outdir/raw/$ltrfilename",
				   -format => 'blast');

    # only one result possible for pairwise bl2seq
    my $result = $ltrin->next_result;
    my $qryname = $result->query_name;                                                               
    my $qrylen  = $result->query_length;
	
    if ($result->hits){
	    
	# only one hit possible
	my $hit = $result->next_hit;
	
	# take only the top hsp
	my $hsp = $hit->next_hsp;
	my $hspeval = $hsp->evalue();                                                             
	my $hspstart = $hsp->start('query');                                                     
	my $hspend = $hsp->end('query');                                                            
	my $hsplen = $hsp->end('query') - $hsp->start('query') + 1;                                     
	my $hspcov = sprintf ("%.2f", ($hsplen / $result->query_length));       
	
	# prune off the LTR based on blast to LTR
	open (TWO, ">/tmp/$id1/$id1.2.rvzfa");
	my $string = $sequence->{'full'}->{$id1};
	my $stringlen = length ($string);
#	my $cleaveme = 20 + $hspend;
#	my $cleaveme = 20 + $ltrlen; # hard coded primer length
#	if ($stringlen >= $cleaveme){
#	    $string =~s/^.{$cleaveme}//;
	if ($stringlen >= $hsplen){
	    $string =~s/^.{$hsplen}//;
	    $sequence->{'noltr'}->{$id1} = $string;
	}
	else {
	    $string = '';
	    $sequence->{'noltr'}->{$id1} = $string;
	}
	print TWO ">$id1\n$string\n";
	close (TWO);

	# short circuit if the pruned sequence is empty
	unless ($string){
	    print LOG "$qryname\t$qrylen\t$hspeval\t$hspstart\t$hspend\t$hsplen\t$hspcov\t";
            print LOG "TOO_SHORT\n";
            print TOOSHORT ">$id1\n$sequence->{'full'}->{$id1}\n";
	    `rm -rf /tmp/$id1/`;
	    next;
	}
		
	# krv blast2seq
	my $krvfilename = blast ($conf, "/tmp/$id1/$id1.2.rvzfa", $refrv, "$outdir/raw");  
	my $krvin = Bio::SearchIO->new(-file   => "$outdir/raw/$krvfilename",
				       -format => 'blast');

	# only one result possible for pairwise bl2seq                                                
	my $kresult = $krvin->next_result;
	my $kqryname = $kresult->query_name;
	my $kqrylen  = $kresult->query_length;

	if ($kresult->hits){
	    
	    # only one hit possible                                                                          
	    my $khit = $kresult->next_hit;

	    # take only the top hsp                                                                      
	    my $khsp = $khit->next_hsp;
	    my $khspeval = $khsp->evalue();
	    my $khspstart = $khsp->start('query');
	    my $khspend = $khsp->end('query');
	    my $khsplen = $khsp->end('query') - $khsp->start('query') + 1;
	    my $khspcov = sprintf ("%.2f", ($khsplen / $kresult->query_length));
	    print LOG "$qryname\t$qrylen\t$hspeval\t$hspstart\t$hspend\t$hsplen\t$hspcov\t";
	    print LOG "$khspeval\t$khspstart\t$khspend\t$khsplen\t$khspcov\n";
	    print RV ">$id1\n$sequence->{'full'}->{$id1}\n";
	    print NOLTRRV ">$id1\n$sequence->{'noltr'}->{$id1}\n";
	}
	else {
	    print LOG "$qryname\t$qrylen\t$hspeval\t$hspstart\t$hspend\t$hsplen\t$hspcov\t";
	    print LOG "NO_RV_HIT\n";
	    print FLANK ">$id1\n$sequence->{'full'}->{$id1}\n";
	    print NOLTRFLANK ">$id1\n$sequence->{'noltr'}->{$id1}\n";
	}
    }
    else {
	print LOG "$qryname\t$qrylen\t";
	print LOG "NO_LTR_HIT\n";
	print NOLTRHITS ">$id1\n$sequence->{'full'}->{$id1}\n";
    }

    # rm the tmp files
    `rm -rf /tmp/$id1`;
    
}
close (RV);
close (FLANK);
close (NOLTRHITS);
close (LOG);

# now find duplicates in the flanking sequences using CD-HIT
`cd-hit-est -i $outdir/flank.fa -o $outdir/flank_cdhit_clust.fa $conf->{'CD-HIT'}`;
`cd-hit-est -i $outdir/flank_noltr.fa -o $outdir/flank_noltr_cdhit_clust.fa $conf->{'CD-HIT'}`;

# mine the flanking sequences
my $insites = {};
my $insites_revcomps = {};
my $insites_revs = {};
my $in = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/flank_noltr_cdhit_clust.fa");
#my $in = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/flank_noltr.fa");
open (I, ">$outdir/insertion.sites");
while (my $sequence_obj = $in->next_seq()){
    my $id  = $sequence_obj->display_id();
    my $seq = $sequence_obj->seq();
    my $len = $sequence_obj->length();

    if ($len < $intlen){
	print I "$id\tTOO_SHORT\n";
	next;
    }
    else {
	my $subseq = substr ($seq, 0, $intlen);
	my $subseqrev = reverse ($subseq);
	my $subseqrevcomp = $subseqrev;
	$subseqrevcomp =~ tr/ACGTacgt/TGCAtgca/;

	print I "$id\t$subseq\t$subseqrev\t$subseqrevcomp\n";
	$insites_revcomps->{$subseq} = $subseqrevcomp;
	$insites->{$subseq}++;
	$insites_revs->{$subseq} = $subseqrev;
    }
}
close (I);

# print the insertion site freqs
open (J, ">$outdir/insertion.sites.freqs");
foreach my $isite (sort keys %$insites){
    print J "$isite\t$insites_revs->{$isite}\t$insites_revcomps->{$isite}\t$insites->{$isite}\n";
}
close (J);

#####SUBS#####

sub blast{
    my $conf   = shift;
    my $onefile  = shift;
    my $twofile  = shift;
    my $out    = shift;

    my $onefile_name;
    if ($onefile =~/\//g){
	$onefile =~m/.*\/(.*)$/;
	$onefile_name = $1;
    }

    else {
	$onefile_name = $onefile;
    }

    my $twofile_name;
    if ($twofile =~/\//g){
        $twofile =~m/.*\/(.*)$/;
        $twofile_name = $1;
    }

    else {
        $twofile_name = $twofile;
    }

    my $cmd = "bl2seq";
    $cmd .= " $conf->{BL2SEQ}";
    $cmd .= " -i $onefile";
    $cmd .= " -j $twofile";
    $cmd .= " -o $out/$onefile_name.$twofile_name.out";
    `$cmd`;
    
    return ("$onefile_name.$twofile_name.out");
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
    close (F);
    return (\%config);
}
