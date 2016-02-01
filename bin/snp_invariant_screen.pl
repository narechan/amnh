#!/usr/bin/perl -w

# This program accepts a 1000G tab snp list
# and screens out invariant characters

# -i is the input snp table

# OR

# -a is the alignment file and -f is the format
# -q if you want to discount missing data and gaps as characters
#    (if missing data is considered invariant)
# -j if you want to discount characters with non-dna residues (N,.. etc)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $alnfile, $format, $dismissing, $nondna);
GetOptions(
    'i|infile=s'  => \$infile,
    'a|aln=s'     => \$alnfile,
    'f|format=s'   => \$format,
    'q|dismissing' => \$dismissing,
    'j|nondna'    => \$nondna,
    );


#####MAIN#####

# readin the file and store data
my $snps = {};

if ($infile){
    open (I, "$infile");
    while (my $line = <I>){
	chomp $line;
	
	my ($acc, $seq) = 
	    split (/\s+/, $line);
	
	# blow apart the sequence
	my @seq = split (//, $seq);
	
	# plug in the snp chars
	my $pos = 0;
	foreach my $snp (@seq){
	    $pos++;
	    $snps->{$pos}->{$acc} = $snp;
	}
#	print STDERR "$acc\t$pos\n";
    }
    close (I);
}

elsif ($alnfile){
    my $alnin = Bio::AlignIO->new(-file   => "$alnfile",
				  -format => "$format");
    my $alnobj = $alnin->next_aln();
    foreach my $seq ($alnobj->each_seq){
        my $id        = $seq->display_id;
        my $sequence = $seq->seq;

	# blow apart the sequence                                                                              
        my @seq = split (//, $sequence);
	
        # plug in the snp chars                                                                              
        my $pos = 0;
        foreach my $snp (@seq){
            $pos++;
            $snps->{$pos}->{$id} = $snp;
        }
#        print STDERR "$id\t$pos\n";
    }
}

else {
    print STDERR "Need alignment data in a legal format!\n";
    die;
}

my $screenedsnps = {};
my $scount = 0;
my $icount = 0;
my $taxa = {};
foreach my $pos (sort {$a <=> $b} keys %$snps){
#    last if ($pos == 67);
    my $uniqbases = {};
    foreach my $acc (sort keys %{$snps->{$pos}}){
	$taxa->{$acc} = 1;

	# screen if missing data is considered invariant
	if ($dismissing){
	    next if ($snps->{$pos}->{$acc} eq "?");
	    next if ($snps->{$pos}->{$acc} eq "-");
	}
	if ($nondna){
	    next unless ($snps->{$pos}->{$acc} =~m/[ACGT]/i);
	}
	
	$uniqbases->{$snps->{$pos}->{$acc}}++;
    }
    my @string;
    my $basecount = 0;
    foreach my $base (sort keys %$uniqbases){
	$basecount++;
	push (@string, $base, $uniqbases->{$base});
    }
    my $string = join "\t", @string;
    print STDERR "$pos\t$basecount\t$string\n";
    
    if ($basecount == 1){
	next;
    }
    else{
	if ($string =~m/[A-Z-]/i){
	    $scount++;
	}
	elsif ($string =~m/[01]/){
	    $icount++;
	}
	else {
	    print STDERR "Unknown Character\n";
#	    die;
	    next;
	}
       
	foreach my $acc (sort keys %{$snps->{$pos}}){
	    push (@{$screenedsnps->{$acc}}, $snps->{$pos}->{$acc});
	}
    }
}

my $taxcount = keys %$taxa;
my $charcount = $scount + $icount;

print " $taxcount $charcount\n";
foreach my $tax (sort keys %$screenedsnps){
    my $seqstring = join "", @{$screenedsnps->{$tax}};
    print "$tax\t$seqstring\n";
}

print STDERR "$scount\t$icount\n";
