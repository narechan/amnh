#!/usr/bin/perl -w

# This program accepts a directory of vcf files
# and creates a matrix from that data for subsequent phylogenetic analysis.

# NOTE: the number of taxa in the matrix is the number of queries in 
#       your list plus 1, for the reference.

# must define  -i and/or -t

# -i is the input vcf directory
# -t is the input table directory (from show-snp in mummer ... prescreen indels)
# -o is the outdir
# -q is the snp quality above which snps are accepted
# -r is the reference (string as in the taxon name)
# -s is a file of ranges to be screened (optional)
# -v is set to on if you want filtered vcfs (optional)
# -p is set to on if you want the number of pairwise snps between all leaves (optional)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($indir, $tabledir, $outdir, $quality, $ref, $screenfile, $vcfson, $pon);
GetOptions(
	   'i|indir=s'  => \$indir,
           't|tabledir=s' => \$tabledir,
	   'o|outdir=s'  => \$outdir,
	   'q|quality=s' => \$quality,
	   'r|reference=s' => \$ref,
           's|screenfile=s' => \$screenfile,
           'v|vcfson' => \$vcfson,
           'p|pon'    => \$pon,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# readin screenfile
my @ints;
if ($screenfile){
    open (L, "$screenfile");
    while (my $line = <L>){
	chomp $line;
	push (@ints, $line);
    }
}


# readin the file and store data
my $snps = {};
my $refbases = {};
my $annots = {};
my @qs;

# make filtered vcfs dir if vcfson is on
(`mkdir -p $outdir/vcfs_filtered`) if ($vcfson);

if ($indir){

    # read in all vcfs and create the query list                                                             
    opendir (D, "$indir");
    my @vcfs = sort (readdir (D));
    shift @vcfs;
    shift @vcfs;
    closedir (D);
    
    #@qs = @vcfs;
    push (@qs, @vcfs);

    foreach my $vcf (@vcfs){
	(open (O, ">$outdir/vcfs_filtered/$vcf")) if ($vcfson);
	
	my $vcfcounter = 0;
	open (I, "$indir/$vcf");
	while (my $line = <I>){
	    chomp $line;
	    next if ($line =~m/^\#/);
	    
	    my ($refctg,
		$refpos,
		$rsid,
		$refbase,
		$snp,
  		$sq,
		$filter,
		$info,
		$format,
		$genotypeinfo) = split (/\t/, $line);
	    my @genotypeinfo;
	    if ($genotypeinfo){
		@genotypeinfo = split (":", $genotypeinfo);
	    }
	    
	    # primary filters
	    next if ($refbase eq "N"); #ambiguous ref                                                           
	    next if ($info =~m/INDEL/); #skip indels 
	    next if ($snp =~m/^I/); #this deals with pacbio vcf files which are version 3.3                   
	    next if ($snp =~m/^D/); #same   

	    # screens
	    my $signal = 0;
	    foreach my $int (@ints){
		my ($start, $end) = split (/\-/, $int);
		if (($refpos >= $start) and ($refpos <= $end)){
		    $signal++;
		}
	    }
	    next if ($signal > 0);
	    
	    # check to see if this is an mnp (vcf 3.3 / pacbio / freebayes handler)
#	    if ((length ($snp) > 1) and ($snp !~m/\,/) and ($info !~m/INDEL/) and ($snp !~m/^I/) and ($snp !~m/^D/)){ #this likely will not work for freebayes!
	    if ((length ($snp) > 1) and ($snp !~m/\,/)){
		my @mnp = split (//, $snp);
		my @ref = split (//, $refbase);
		my $mnpcount = @mnp;
		for (my $i = 0; $i < $mnpcount; $i++){
		    
		    # plug in the reference chars                                                           
		    $snps->{$refctg}->{$refpos}->{$ref}  = $ref[$i];
		    $refbases->{$refctg}->{$refpos} = $ref[$i];

		    # secondary filters                                                                  
		    if ($sq < $quality){ #sq                                                            
			$snps->{$refctg}->{$refpos}->{$vcf} = "?";
			next;
		    }
		    if ($snp =~m/\,/){ #call hets                                                         
			$snps->{$refctg}->{$refpos}->{$vcf} = "?";
			next;
		    }
		    if ($genotypeinfo){
			if ($genotypeinfo[0] eq "0/1"){ #genotype hets                                        
			    $snps->{$refctg}->{$refpos}->{$vcf} = "?";
			    next;
			}
		    }

		    # plug in the snp chars                                                             
		    $snps->{$refctg}->{$refpos}->{$vcf} = $mnp[$i];
		    $vcfcounter++;
		    $refpos++;
		    
		    # note that for mnps, the vcf line will be replicated
		    # in accordance with the length of the mnp
		    (print O "$line\n") if ($vcfson);
		}
	    }
	    else{
		# plug in the reference chars
		$snps->{$refctg}->{$refpos}->{$ref}  = $refbase;
		$refbases->{$refctg}->{$refpos} = $refbase;
		
		# secondary filters
		if ($sq < $quality){ #sq
		    $snps->{$refctg}->{$refpos}->{$vcf} = "?";
		    next;
		}
		if ($snp =~m/\,/){ #call hets
		    $snps->{$refctg}->{$refpos}->{$vcf} = "?";
		    next;
		}
		if ($genotypeinfo){
		    if ($genotypeinfo[0] eq "0/1"){ #genotype hets
			$snps->{$refctg}->{$refpos}->{$vcf} = "?";
			next;
		    }
		}

		# plug in the snp chars
		$snps->{$refctg}->{$refpos}->{$vcf} = $snp;
#		print STDERR "$refpos\n";
		$vcfcounter++;
		(print O "$line\n") if ($vcfson);
	    }
	}
	print STDERR "$vcf\t$vcfcounter\n";
	close (I);
	(close (O)) if ($vcfson);
    }
}
if ($tabledir){
    
    # read in all tables and create the query list                                                         
    opendir (D, "$tabledir");
    my @tables = sort (readdir (D));
    shift @tables;
    shift @tables;
    closedir (D);

#    @qs = @tables;
    push (@qs, @tables);
    
    foreach my $table (@tables){

        my $tablecounter = 0;
        open (I, "$tabledir/$table");
        while (my $line = <I>){
            chomp $line;
            next if ($line =~m/^\#/);

	    my ($refctg, 
		$query, 
		$refpos, 
		$refbase, 
		$snp, 
		$sq) = split (/\t/, $line);

            # primary filters (note that indels should be pre-screened)
            next if ($refbase eq "N"); #ambiguous ref                                                          

	    # screens                                                                                         
            my $signal = 0;
            foreach my $int (@ints){
                my ($start, $end) = split (/\-/, $int);
                if (($refpos >= $start) and ($refpos <= $end)){
                    $signal++;
                }
            }
            next if ($signal > 0);

            # plug in the reference chars                                                                  
            $snps->{$refctg}->{$refpos}->{$ref}  = $refbase;
            $refbases->{$refctg}->{$refpos} = $refbase;

            # secondary filters                                                                              
            if ($sq < $quality){ #sq                                                                         
                $snps->{$refctg}->{$refpos}->{$table} = "?";
                next;
            }
            if ($snp =~m/\,/){ #hets --> only way to do in table format
                $snps->{$refctg}->{$refpos}->{$table} = "?";
                next;
            }
	    # plug in the snp chars                                                                          
            $snps->{$refctg}->{$refpos}->{$table} = $snp;
            $tablecounter++;
        }
        print STDERR "$table\t$tablecounter\n";
        close (I);
    }
}



# post process datastruc to fill in ref alleles
# for queries with snp gaps
foreach my $ctg (keys %$snps){
    foreach my $pos (keys %{$snps->{$ctg}}){
	foreach my $query (@qs){
	    if (exists($snps->{$ctg}->{$pos}->{$query})){
		next;
	    }
	    else {
		$snps->{$ctg}->{$pos}->{$query} = $refbases->{$ctg}->{$pos};
	    }
	}
    }
}

# get the pairwise counts of diffs between queries
# and print to stderr
if ($pon){
    my $pairwise = {};
    my $qrycopy = {};
    my $qrycopy2 = {};
    foreach my $qry (@qs){
	$qrycopy->{$qry} = 1;
	$qrycopy2->{$qry} = 1;
    }
    
    foreach my $q1 (sort keys %$qrycopy){
	foreach my $q2 (sort keys %$qrycopy){
	    next if ($q1 eq $q2);
	    
#	print STDERR "$q1:$q2\n";
	    my $counter = 0;
	    my $pair = $q1 . "-" . $q2;
	    foreach my $ctg (sort keys %$snps){
		foreach my $pos (sort {$a <=> $b} keys %{$snps->{$ctg}}){
#		print STDERR "$counter\t$ctg\t$pos\n";
		    if ( ($snps->{$ctg}->{$pos}->{$q1} eq "?") or ($snps->{$ctg}->{$pos}->{$q2} eq "?") ){
			next;
		    }
		    elsif ($snps->{$ctg}->{$pos}->{$q1} ne $snps->{$ctg}->{$pos}->{$q2}){
			$pairwise->{$pair}++;
			$counter++;
#		    print STDERR "$counter\t$ctg\t$pos\n";
		    }
		    else {
			next;
		    }
		}
	    }
	}
	delete $qrycopy->{$q1};
    }
    
    foreach my $q1 (sort keys %$qrycopy2){
	foreach my $q2 (sort keys %$qrycopy2){
	    next if ($q1 eq $q2);
	    my $pair = $q1 . "-" . $q2;
	    if (exists ($pairwise->{$pair})){
		print STDERR "$pair\t$pairwise->{$pair}\n";
	    }
	    else {
		print STDERR "$pair\t0\n";
	    }
	}
	delete $qrycopy2->{$q1};
    }
}

# get some matrix data
my $taxa = @qs;
$taxa++;
my $alnlentot = 0;
foreach my $c (keys %$snps){
    foreach my $s (keys %{$snps->{$c}}){
	$alnlentot++;
    }
}

# sort print the matrix and charpars                                   
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
open (TBL, ">$outdir/table");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

# print table header
foreach my $contig (keys %$snps){
    print TBL "Contig\tSNP\t";
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	my @sp;
	foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    push (@sp, $sp);
	}
	my $spstring = join "\t", @sp;
	print TBL "$spstring\n";
	last;
    }
    last;
}

my $start = 1;
my $end;
#my $vreseen = 0; ##
#my $vrenew  = 0; ##
foreach my $contig (keys %$snps){
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	$end = $start;
	print CAT "[Partition ctg$contig snp$count]\n";

	my @table;
	push (@table, $contig, $count);
	
	my $charset = "ctg" . $contig . "snp" . $count;
	print PRT "CHARSET $charset=$start-$end;\n";
    
	my @qrysnps;
	foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    print CAT "$sp\t$snps->{$contig}->{$count}->{$sp}\n"; 
	    push (@table, $snps->{$contig}->{$count}->{$sp});
	    push (@qrysnps, $snps->{$contig}->{$count}->{$sp}) unless ($sp eq $ref);
	}

=head
	## vre specific hack ##
	my @vreqrysnps;
	my $vresignal = 0;
	my $poop = $snps->{$contig}->{$count}->{'VRE1107.vcf'};
        foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    next if ($sp eq "VRE1107.vcf");
	    if ($poop eq $snps->{$contig}->{$count}->{$sp}){
		$vresignal++;
	    }
	}
	if ($vresignal > 0){
	    $vreseen++;
	}
	else{
	    $vrenew++;
	}

	## end vre specefic hack ##
=cut
	
	# check to see if all the qry snps are the same
	if (keys %{{ map {$_, 1} @qrysnps }} == 1){
	    push (@table, "N");
	}
	else {
	    my $signal = 0;
	    foreach my $qsnp (@qrysnps){
		($signal++) if ($qsnp eq "?");
	    }
	    if ($signal > 0){
		push (@table, "*");
	    }
	    else {
		push (@table, "+");
	    }
	}
	my $table = join "\t", @table;
	print TBL "$table\n";
#	print TBL "$table\t$annots->{$contig}->{$count}\n";

	print CAT "\n";
	$start = $end + 1;
    }
}
print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";
close (CAT);
close (PRT);
close (TBL);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
#print STDERR "$vrenew\t$vreseen\n";
