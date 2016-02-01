#!/usr/bin/perl -w

# This program accepts a directory of freebayes vcf files
# and creates a table of output

# -i is the input vcf directory
# -q is the snp quality above which a position's calls are accepted
# -r is the reference (string as in the taxon name)

#####SETUP#####

use strict;
use Getopt::Long;
use Math::Round;

my ($indir, $quality, $ref);
GetOptions(
	   'i|indir=s'  => \$indir,
	   'q|quality=s' => \$quality,
	   'r|reference=s' => \$ref,
	   );

#####MAIN#####

# read in all vcfs and create the query list
opendir (D, "$indir");
my @vcfs = sort (readdir (D));
shift @vcfs;
shift @vcfs;
closedir (D);

# readin the file and store data
my $snps = {};
my $quals = {};
my $refbases = {};

foreach my $vcf (@vcfs){

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

	# primary filters
	next if ($sq < $quality);
	
	$vcfcounter++;

	# parse info
	my $type;
	my $qa;
	my $ro;
	my $dp;
	my @info = split (/\;/, $info);
	foreach my $i (@info){
	    my ($key, $value) = split (/\=/, $i);
	    if ($key eq "TYPE"){
		$type = $value;
	    }
	    elsif ($key eq "QA"){
		$qa = $value;
	    }
	    elsif ($key eq "DP"){
		$dp = $value;
	    }
	    elsif ($key eq "RO"){
		$ro = $value;
	    }
	    else{
		next;
	    }
	}
	my @type = split (/\,/, $type);
	my @qa   = split (/\,/, $qa);
	
	my @qual;
	foreach my $q (@qa){
	    my $quote = nearest (0.001, $q / $dp);
	    push (@qual, $quote);
	}
	
	# secondary filters
	if ($snp =~m/\,/){ #hets
	    my @snps = split (/\,/, $snp);
	    my $scount = 0;
	    foreach my $s (@snps){
		$scount++;
		next unless ($type[$scount-1] eq "snp"); #skip unless this call is a snp

		$snps->{$refctg}->{$refpos}->{$vcf}->{$scount} = $s;
		$quals->{$refctg}->{$refpos}->{$vcf}->{$scount} = $qual[$scount-1];

		# plug in the reference chars                                                          
		$snps->{$refctg}->{$refpos}->{$ref}->{1}  = $refbase;
#		$quals->{$refctg}->{$refpos}->{$ref}->{1} = nearest (0.001, $ro / $dp);
		$refbases->{$refctg}->{$refpos} = $refbase;
		print "$line\n";
	    }
	}
	else {
	    next unless ($type eq "snp"); #skip unless we're dealing with a snp
	    $snps->{$refctg}->{$refpos}->{$vcf}->{1} = $snp;
	    $quals->{$refctg}->{$refpos}->{$vcf}->{1} = $qual[0];
	    
	    # plug in the reference chars                                                                   
	    $snps->{$refctg}->{$refpos}->{$ref}->{1}  = $refbase;
#	    $quals->{$refctg}->{$refpos}->{$ref}->{1} = nearest (0.001, $ro / $dp);
	    $refbases->{$refctg}->{$refpos} = $refbase;
	    print "$line\n";
	}
    }
    print STDERR "$vcf\t$vcfcounter\n";
    close (I);
}

# print out report with a line per position
my $pooplocal = 0;
my $culturelocal = 0;
my $mixed = 0;
foreach my $contig (keys %$snps){
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){

        my @table;
        push (@table, $contig, $count);

	# vre specefic data
	my $poop = $snps->{$contig}->{$count}->{'VRE1107_ploidy10_pooledBothWays_F0_C2.vcf'};
        my $efdo = $snps->{$contig}->{$count}->{'EFDO'};
	
	my $spcount = 0;
	foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
	    next if ($sp eq "VRE1107_ploidy10_pooledBothWays_F0_C2.vcf");
	    next if ($sp eq "EFDO");
	    if (exists ($snps->{$contig}->{$count}->{$sp})){
		$spcount++;
	    }
	    else{
		next;
	    }
	}
	
	if (($poop) and $spcount == 0){
	    $pooplocal++;
	}
	elsif (( !defined($poop)) and ($spcount > 0)){
	    $culturelocal++;
	}
	else {
	    $mixed++;
	}
    }
}

print STDERR "$pooplocal\t$culturelocal\t$mixed\n";
	    
# post process datastruc to fill in ref alleles                                                            
# for queries with snp gaps                                                                             
foreach my $ctg (keys %$snps){                                                                          
    foreach my $pos (keys %{$snps->{$ctg}}){                                                                
        foreach my $query (@vcfs){                                                                          
            if (exists($snps->{$ctg}->{$pos}->{$query}->{1})){                                               
                next;                                                                                     
            }                                                                                               
            else {                                                                                        
                $snps->{$ctg}->{$pos}->{$query}->{1} = $refbases->{$ctg}->{$pos};                           
            }                                                                                               
        }                                                                                                   
    }                                                                                                    
}                       	

foreach my $contig (keys %$snps){
    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	
        my @table;
        push (@table, $contig, $count);
        
	my @qrysnps;
        foreach my $sp (sort keys %{ $snps->{$contig}->{$count} }){
#	    print STDERR "$sp\n";
	    my @all;
	    foreach my $allele (sort keys %{ $snps->{$contig}->{$count}->{$sp} }){
		my $a = $snps->{$contig}->{$count}->{$sp}->{$allele};
		my $b = $quals->{$contig}->{$count}->{$sp}->{$allele};
		my $string;
		if ($b){
		    $string = $a . "(" . $b . ")";
		}
		else{
		    $string = $a . "(" . "NA" . ")";
		}
		push (@all, $string);
	    }
	    my $all = join ";", @all;
            push (@table, $all);
            push (@qrysnps, $all) unless (($sp eq $ref) or ($sp eq "VRE1107_ploidy10_pooledBothWays_F0_C2.vcf"));
        }
	my $table = join "\t", @table;
#        print "$table\n";
	
    }
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
    }
}
=cut

=head
# post process datastruc to fill in ref alleles
# for queries with snp gaps
foreach my $ctg (keys %$snps){
    foreach my $pos (keys %{$snps->{$ctg}}){
	foreach my $query (@vcfs){
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
my $pairwise = {};
my $qrycopy = {};
my $qrycopy2 = {};
foreach my $qry (@vcfs){
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


# get some matrix data
my $taxa = @vcfs;
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
my $vreseen = 0; ##
my $vrenew  = 0; ##
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
print STDERR "$vrenew\t$vreseen\n";
=cut
