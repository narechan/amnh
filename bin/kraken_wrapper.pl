#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('a:b:f:d:o:t:', \%opts);
my $afile  = $opts{'a'}; #fastq1
my $bfile  = $opts{'b'}; #fastq2
my $ffile  = $opts{'f'};  #concat fasta if you already have it
my $db     = $opts{'d'}; #kraken db
my $tlookup = $opts{'t'}; #taxonomy lookup
my $outdir  = $opts{'o'};

`mkdir -p $outdir`;

# either a/b or f need to be specified

# concatenate reads and add ambiguous character so no artificial
# kmers are processed.
if ($afile and $bfile){
    my $cat = {};
    
    my $acc;
    my $desc;
    my $seq;
    my $acounter = 0;
    open (A, "$afile");
    while (my $line = <A>){
	chomp $line;
	$acounter++;
	
	if ($acounter == 1){
	    ($acc, $desc) = split (/\s/, $line);
	    $acc =~s/\@//g;
	}
	elsif ($acounter == 2){
	    $seq = $line;
	}
	elsif ($acounter == 4){
	    $cat->{$acc} = $seq;
	    $acounter = 0;
	}
	else {
	    next;
	}
    }
    
    my $acc;
    my $desc;
    my $seq;
    my $bcounter = 0;	
    open (B, "$bfile");
    while (my $line = <B>){
	chomp $line;
	$bcounter++;
	
	if ($bcounter == 1){
	    ($acc, $desc) = split (/\s/, $line);
	    $acc =~s/\@//g;
	}
	elsif ($bcounter == 2){
	    $seq = $line;
	}
	elsif ($bcounter == 4){
	    $cat->{$acc} .= "N$seq";
	    $bcounter = 0;
	}
	else {
	    next;
	}
    }
    
    open (F, ">$outdir/fasta");
    foreach my $acc (keys %$cat){                                                                            
	print F ">$acc\n$cat->{$acc}\n";                                                                      
    }   
    close (F);

    # run kraken                                                                                           
    `kraken --threads 40 -db "$db" $outdir/fasta > $outdir/kraken.classified`;

}
elsif ($ffile) {
    `kraken --threads 40 -db "$db" $ffile > $outdir/kraken.classified`;
}

else {
    print STDERR "Need some sequence!\n";
    die;
}

# map the taxonomy names to the kraken output file
my $tax = {};
open (T, "$tlookup");
while (my $line = <T>){
    chomp $line;
    $line =~s/\t//g;
    my @line = split (/\|/, $line);
    if ($line[3] eq "scientific name"){
        $tax->{$line[0]} = $line[1];
    }
    else{
        next;
    }
}
close (T);

# parse and add to the kraken output
open (K, "$outdir/kraken.classified");
open (O, ">$outdir/kraken.classified.tname");
while (my $line = <K>){
    chomp $line;
    my @line = split (/\t/, $line);
    my $t = $tax->{$line[2]};
    print O "$line\t$t\n";
}
close (K);
close (O);



=head                                                                                                             
    my $cat = {};                                                                                                     
my $aseqobj = Bio::SeqIO->new(-file   =>"$afile",                                                                 
                              -format =>"fastq");                                                                 
my $counter = 0;                                                                                                  
while (my $aseq = $aseqobj->next_seq()){                                                                          
    $counter++;                                                                                                   
#    print STDERR "$counter\n";                                                                                   
    my $aid            = $aseq->display_id();                                                                     
    my $asequence      = $aseq->seq();                                                                            
    $cat->{$aid} = $asequence;                                                                                    
}                                                                                                                 
                                                                                                                  
my $bseqobj = Bio::SeqIO->new(-file   =>"$bfile",                                                                 
                              -format =>"fastq");                                                                 
while (my $bseq = $bseqobj->next_seq()){                                                                          
    my $bid            = $bseq->display_id();                                                                     
    my $bsequence      = $bseq->seq();                                                                            
    $cat->{$bid} .= "N$asequence";                                                                                
}                                                                                                                 
                                                                                                                  
foreach my $acc (keys %$cat){                                                                                     
    print ">$acc\n$cat->{$acc}\n";                                                                                
}                                                                                                                 
=cut                                                                                                              

