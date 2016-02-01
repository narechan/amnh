#!/usr/bin/perl -w

=head1 NAME

txt_slurp_aln_rand.pl

=head1 SYNOPSIS

  txt_slurp_aln_rand.pl -- this program will partition a clean text
    string  into "gene" sizes of arbitrary length (select the first x characters for the first gene, 
    first taxa; select the second  x characters for the first gene, second taxa; etc...) and align them, 
    or accept an alignment and randomize in 1 of 2 ways:

    1. randomize rows (scramble each taxon's alignment) 
    2. randomize columns (scramble taxa foreach position)

Options:

 --help        Show brief help and exit
 --infile      Contains the text to be analyzed (single string slurp) (optional)
 --alignfile   Contains the alignment to be randomized (optional)
 --config      Is the configuration for your run (MAFFT params)
 --chars       The number of chars in each "gene"
 --total       The total number of chars you want
 --members     The number of taxa members to simulate randomly
 --outdir      Is your output dir
 --method      Is the randomization method you want to employ of the 2 listed above (1 or 2)
    if not specified or specified with some other value, then it outputs the unrandomized matrix (optional)
 --procs       Is the number of parallel processes you want to run
 --scramble    Is set if you want select partitions from the whole in a scrambled way; if not 
    set it will select partitions of size specified in --chars from the top down as the infile or
    alignment is read. (optional)

mafft must be in your path

The config must specify parameters for the mafft run

=head1 DESCRIPTION

Run mafft on extracted "genes" from a "genomic" text, and random taxa

Usage examp:

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
use Bio::SeqIO;
use Bio::AlignIO;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Parallel::ForkManager;

my ($help, $infile, $alnfile, $config, $outdir, $chars, $total, $members, $method, $procs, $scramble);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'a|alignfile=s'   => \$alnfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'x|chars=s'       => \$chars,
    't|total=s'       => \$total,
    'm|members=s'     => \$members,
    'y|method=s'      => \$method,
    'p|procs=s'       => \$procs,
    's|scramble'      => \$scramble,
    ) or pod2usage;

pod2usage if $help;

for my $option ($outdir, $config, $chars, $total, $members, $method, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/alns`;
`mkdir -p $outdir/subfastas`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

my $data = {};

# txt slurp
if ($infile){
    my $txt;
    open (F, "$infile");
    while (my $line = <F>){
	chomp $line;
	$line =~s/\s//g;
	$line =~s/\d//g;
	$line =~tr/[a-z]/[A-Z]/;
	$line =~s/[[:punct:]]//g;
	$txt .= $line;
    }
    close (F);
    
    # if splitting out chunks and aligning subsections
    # generate randome subfastas of requested size over the entire text
    my $counter = 0;
    while (1){
        $counter++;
	
	my $textleftlen = length $txt;
	last if ($textleftlen < ($members * $chars));
	
	open (O, ">$outdir/subfastas/part$counter");
        for (my $i = 1; $i <= $members; $i++){
            my $substring = substr ($txt, 0, $chars, '');
            print O ">tax$i\n$substring\n";
        }
        close (O);

    }
    
    # align all subfastas
    opendir (A, "$outdir/subfastas");
    my @fastas = sort readdir (A);
    shift @fastas;
    shift @fastas;
    closedir (A);
    my $pm = Parallel::ForkManager->new($procs);
    foreach my $fasta (sort @fastas){
	$pm->start and next;

	print STDERR "Mafft $fasta\n";
	mafft_run ($fasta, $conf, $outdir);

	$pm->finish;
    }
    $pm->wait_all_children;
    
    # store and concatenate the aln data
#    opendir (J, "$outdir/alns");
#    my @alns = sort readdir (J);
#    shift @alns;
#    shift @alns;
#    closedir (J);

    my $concataln = {};
#    foreach my $aln (@alns){
    for (my $i = 1; $i <= $counter - 1; $i++){
	my $alnin = Bio::AlignIO->new(-file   => "$outdir/alns/part$i",
                                      -format => "fasta");
	my $alnobj = $alnin->next_aln();
        foreach my $seq ($alnobj->each_seq){
            my $id        = $seq->display_id;
            my $sequence = $seq->seq;
	    
	    $concataln->{$id} .= $sequence;
	}
    }

    open (B, ">$outdir/alns/bigconcat.aln");
    foreach my $acc (sort keys %$concataln){
	print B ">$acc\n$concataln->{$acc}\n";
    }
    close (B);
    
    # break apart the alignment and (randomly) select partitions
    my $alnin = Bio::AlignIO->new(-file   => "$outdir/alns/bigconcat.aln",
				  -format => "fasta");
    my $alnobj = $alnin->next_aln();
    my $alnlen = $alnobj->length;

    my $start = 1;
    my $alndata = {};
    my $tcounter = 0;
    my $taxa = {};
    my @partnames;
    for (my $start = 1; $start < $alnlen; $start = $start + $chars){
        $tcounter++;

        my $end;
        if (($start + $chars - 1) > $alnlen){
#            $end = $alnlen;
	    last; #get rid of the last truncated partition
        }
        else{
            $end = $start + $chars - 1;
        }

        my $partname = 'part' . $tcounter;
        push (@partnames, $partname);

        foreach my $seq ($alnobj->each_seq){
            my $id        = $seq->display_id;
            my $partition = $seq->subseq($start, $end);

            $taxa->{$id} = 1;
            $alndata->{$partname}->{$id} = $partition;
        }
    }
    
    # randomly (or not) select partitions from the ones generated to create the final matrix         
    my $partsneeded = $total / $chars;
    if ($scramble){
	for (my $j = 1; $j <= $partsneeded; $j++){
	    my $p = $partnames[int(rand(@partnames))];
	    if (exists ($data->{$p})){
		($p = $partnames[int(rand(@partnames))]) until (! exists($data->{$p}));
	    }
	    foreach my $id (sort keys %{$alndata->{$p}}){
		my $pa = $alndata->{$p}->{$id};
		my @pa = split (//, $pa);
		$data->{$p}->{$id} = [@pa];
	    }
	}
    }
    else {
#	foreach my $p (@partnames){
	for (my $j = 1; $j <= $partsneeded; $j++){
	    my $p = "part" . $j;
	    foreach my $id (sort keys %{$alndata->{$p}}){
                my $pa = $alndata->{$p}->{$id};
                my @pa = split (//, $pa);
                $data->{$p}->{$id} = [@pa];
            }
        }
    }
    print yellow;
    
}

# real aln: create artificial partitions from the concatenated alignment
# randomly select from the total alignment the subalignments of specified length
# no dups allowed
elsif ($alnfile){
    my $alnin = Bio::AlignIO->new(-file   => "$alnfile",
				  -format => "nexus");
    my $alnobj = $alnin->next_aln();
    my $alnlen = $alnobj->length;
    
    my $start = 1;
    my $alndata = {};
    my $counter = 0;
    my $taxa = {};
    my @partnames;
    for (my $start = 1; $start < $alnlen; $start = $start + $chars){
	$counter++;

	my $end;
	if (($start + $chars - 1) > $alnlen){
#	    $end = $alnlen;
	    last; #get rid of the last truncated partition
	}
	else{
	    $end = $start + $chars - 1;
	}

	my $partname = 'part' . $counter;
	push (@partnames, $partname);

	foreach my $seq ($alnobj->each_seq){
	    my $id        = $seq->display_id;
	    my $partition = $seq->subseq($start, $end);

	    $taxa->{$id} = 1;
	    $alndata->{$partname}->{$id} = $partition;
	}
    }

    # get a random unique subset of the taxa for this matrix                          
    my @taxa = keys %$taxa;
    my $species = {};
    for (my $j = 1; $j <= $members; $j++){
        my $tax = $taxa[int(rand(@taxa))];
        if (exists ($species->{$tax})){
            ($tax = $taxa[int(rand(@taxa))]) until (! exists($species->{$tax}));
        }
        $species->{$tax} = 1;
    }
    
    # randomly (or not) select partitions from the ones generated to create the final matrix
    my $partsneeded = $total / $chars;
    if ($scramble){
	for (my $j = 1; $j <= $partsneeded; $j++){
	    my $p = $partnames[int(rand(@partnames))];
	    if (exists ($data->{$p})){
		($p = $partnames[int(rand(@partnames))]) until (! exists($data->{$p}));
	    }
	    
	    my $tc = 0;
	    foreach my $id (sort keys %$species){
		$tc++;
		my $pa = $alndata->{$p}->{$id};
		my @pa = split (//, $pa);
#		$data->{$p}->{$id} = [@pa];
		my $newid = "tax". $tc;
		$data->{$p}->{$newid} = [@pa]
	    }
	}
    }
    else{
#	foreach my $p (@partnames){
	for (my $j = 1; $j <= $partsneeded; $j++){
            my $p = "part" . $j;
	    my $tc = 0;
            foreach my $id (sort keys %$species){
                $tc++;
                my $pa = $alndata->{$p}->{$id};
                my @pa = split (//, $pa);
#               $data->{$p}->{$id} = [@pa];                                                                 
                my $newid = "tax". $tc;
                $data->{$p}->{$newid} = [@pa]
            }
	}
    }
    print yellow;
}
else {
    print STDERR "Need an allowable input file\n";
    die;
}

# Build the matrices after randomization method, if any.
my $newdata = {};
print STDERR "Building matrix\n";
if ($method == 1){
    foreach my $part (sort keys %$data){
	my $taxcounter = 0;
	foreach my $tax (sort keys %{$data->{$part}}){
#	    $taxcounter++;
#	    if ($taxcounter == 1){
#		my @seq = @{$data->{$part}->{$tax}};
#		my $seqstr = join "", @seq;
#		$newdata->{$part}->{$tax} = $seqstr
#	    }
#	    else {
	    my @seq = @{$data->{$part}->{$tax}};
	    my @shuffseq = shuffle @seq;
	    my $shuffseqstr = join "", @shuffseq;
	    $newdata->{$part}->{$tax} = $shuffseqstr;
#	    }
	}
    }
}
elsif ($method == 2){
    # take vertical slices of the matrix, randomize those columns and                                         
    # affix back onto a new, column randomized matrix                                                        
    foreach my $gene (sort keys %$data){
        for (my $j = 0; $j <= $chars - 1; $j++){
            my @column;
	    my @tid;
            foreach my $tid (sort keys %{$data->{$gene}}){
                push (@column, $data->{$gene}->{$tid}[$j]);
		push (@tid, $tid);
            }

            my @columnshuffled = shuffle (@column);
            my $ccounter = -1;
	    foreach my $t (@tid){
		$ccounter++;
#		push @{$newdata->{$gene}->{$t}}, $columnshuffled[$ccounter];
#		print STDERR "$gene\t$t\n";
		$newdata->{$gene}->{$t} .= $columnshuffled[$ccounter];
            }
        }
    }
}
else {
    foreach my $part (sort keys %$data){
	foreach my $tax (sort keys %{$data->{$part}}){
	    my @seq = @{$data->{$part}->{$tax}};
	    my $seq = join "", @seq;
	    $newdata->{$part}->{$tax} = $seq;
	}
    }
}
    
# sort print the concatenation and charpars
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$members NCHAR=$total;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIJKLMNOPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort keys %$newdata){
    $end = $start - 1 + $chars;
    print CAT "[Partition $count length $chars chars $start-$end]\n";
    print PRT "CHARSET $count=$start-$end;\n";
    foreach my $sp (sort keys %{ $newdata->{$count} }){
       print CAT "$sp\t$newdata->{$count}->{$sp}\n";                                                       
    }
    print CAT "\n";
    $start = $end + 1;
}

print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";

close (CAT);
close (PRT);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;

#####SUBS#####

sub mafft_run{
    my $fas  = shift;
    my $conf = shift;
    my $out  = shift;

    # mafft
    my $mafft = "mafft";
    ($mafft .= " $conf->{'MAFFT'}") if ($conf->{'MAFFT'});
    $mafft .= " $out/subfastas/$fas";

    `$mafft 1>$out/alns/$fas 2>$out/stderr/$fas.stderr`;
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

