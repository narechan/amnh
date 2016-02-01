#!/usr/bin/perl -w

# -o is the outdir
# -s is the directory for SNVs (turns snps on)
# -i if you want indels (on/off switch)
# -c is the directory for cnvs (turns cnvs on)

#####SETUP#####

use strict;
use Getopt::Long;
use Text::CSV;

my ($outdir, $snpsdir, $indelsyes, $cnvsdir);
GetOptions(
    'o|outdir=s'  => \$outdir,
    's|snpsdir=s' => \$snpsdir,
    'i|indelsyes' => \$indelsyes,
    'c|cnvsdir=s' => \$cnvsdir,
    );

`mkdir -p $outdir`;

#####MAIN#####

# GLOBALS
my $snps = {};
my $indels = {};
my $pindels = {};
my $refbases = {};
my $class = {};
my $cnvs = {};
my $pcnvs;
my $cnvclass = {};
my $qsnvs = {};
my $qcnvs = {};

if ($snpsdir){
    
    # read in the snp files
    # and store the data
    opendir (D, "$snpsdir");
    my @sfiles = sort (readdir (D));
    shift @sfiles;
    shift @sfiles;
    closedir (D);
    
    foreach my $sfile (@sfiles){
	print STDERR "Processing $sfile\n";

	$qsnvs->{$sfile} = 1;
	my $path = $snpsdir . "/" . $sfile;
    
	open (I, "$path");
	my $cntins = 0;
	my $cntdel = 0;
	my $cntsnp = 0;
	my $cntmul = 0;
	while (my $line = <I>){
	    chomp $line;

            # danny new snp method without blood and more stringent                                         
            next if ($line =~m/^genesymbol/);                                                                   
            $line =~s/\"//g;                                                                                
            my @line = split (/\t/, $line);                                                              
            my $refctg  = $line[0];                                                                     
            my $start   = $line[1];                                                                        
            my $end     = $line[2];                                                                        
            my $refbase = $line[3];                                                                          
            my $snp     = $line[4];                                                                      
            my $refbaselen = length ($refbase);                                                            
            my $snplen     = length ($snp);            

=head
	    # danny new snp method without blood and more stringent
	    next if ($line =~m/^\"start/);
	    $line =~s/\"//g;
	    my @line = split (/\t/, $line);
	    my $refctg  = $line[5];                                                                       
            my $start   = $line[0];                                                                       
            my $end     = $line[6];                                                                   
            my $refbase = $line[7];                                                                   
            my $snp     = $line[8];                                                                        
            my $refbaselen = length ($refbase);                                                         
            my $snplen     = length ($snp);      
=cut
=head
            # danny raw files
	    next if ($line =~m/^Func\tGene/);
	    next if ($line =~m/^\t/);

	    # use TXT::csv to remove tabs between quoted strings
	    my $csv = Text::CSV->new({eol => "\n", sep_char => "\t", binary => 1});
	    $csv->parse($line) or die "csv error: " . $csv->error_input();
	    my @columns = $csv->fields();    
	    s/\t//g for (@columns);
	    $csv->combine(@columns) or die "csv error: " . $csv->error_input();
	    my $string = $csv->string();
	    
	    my @line    = split (/\t/, $string);
	    my $func    = $line[0];
	    my $snpstring = $line[3];
	    my $refctg  = $line[8];
	    my $start   = $line[9];
	    my $end     = $line[10];
	    my $refbase = $line[11];
	    my $snp     = $line[12];
	    my $refbaselen = length ($refbase);
	    my $snplen     = length ($snp);
	    
	    # filters
#	    next unless ($func =~m/exonic/);
#	    next if ($func =~m/ncRNA_exonic/);
=cut
=head
	    # yale raw files
	    next if ($line =~m/^chromosome/);
	    next if ($line =~m/^\n/);
	    next if ($line =~m/^\t/);
	    my @line    = split (/\t/, $line);
	    my $func   = $line[7];
	    my $refctg = $line[0];
	    my $start  = $line[2];
	    my $r      = $line[4];
	    $r =~s/.*\///g;
	    my $refbase = $r;
	    my $s      = $line[5];
	    $s =~s/.*\///g;
	    my $snp = $s;
	    my $refbaselen = length ($r);                                                             
            my $snplen     = length ($s); 
	    
	    my $adder = 0;
	    if ($refbaselen >= $snplen){
		$adder = $refbaselen - 1;
	    }
	    else {
		$adder = $snplen - 1;
	    }
	    my $end = $start + $adder;

#	    print STDERR "$line\n";
	    
	    # screen for non-dna characters and log
	    unless (($refbase=~m/\A[ACGT-]+\z/i) and ($snp =~m/\A[ACGT-]+\z/i)){
		print STDERR "Illegal characters: $line\n";
		next;
	    }
=cut	
	    # common to both file formats from here down
	    # deal with inserts
	    if ($refbase eq "-"){
		$cntins++;
		my $int = $start . "-" . $end;
		$class->{$refctg}->{$int} = "INS";
		$indels->{$refctg}->{$int}->{$sfile} = $snp;
	    }
	    
	    # deal with deletions
	    elsif ($snp eq "-"){
		$cntdel++;
		my $int = $start . "-" . $end;
		$class->{$refctg}->{$int} = "DEL"; 
		$indels->{$refctg}->{$int}->{$sfile} = $refbase;
	    }
	    
	    # deal with multibase changes
#	    elsif (($refbaselen > 1) and ($snplen > 1)){
	    elsif (($refbaselen > 1) or ($snplen > 1)){
		$cntmul++;
		my $int = $start . "-" . $end;
		$class->{$refctg}->{$int} = "MUL";
		$indels->{$refctg}->{$int}->{$sfile} = $snp;
	    }
	    
	    # deal with snps
	    else {
		$cntsnp++;
		$class->{$refctg}->{$start} = "SNP";
		
		# store the reference chars
#		$snps->{$refctg}->{$start}->{$ref}  = $refbase;
		$refbases->{$refctg}->{$start} = $refbase;
		
		# plug in the snp chars
		$snps->{$refctg}->{$start}->{$sfile} = $snp;
	    }
	}
	close (I);
	print STDERR "$cntins\t$cntdel\t$cntmul\t$cntsnp\n";
    }
}

if ($cnvsdir){

    # read in the cnv files
    # and store the data                                                                                     
    opendir (C, "$cnvsdir");
    my @cfiles = sort (readdir (C));
    shift @cfiles;
    shift @cfiles;
    closedir (C);

    foreach my $cfile (@cfiles){
        print STDERR "Processing $cfile\n";
    
	open (I, "$cnvsdir/$cfile");
        $qcnvs->{$cfile} = 1;

        my $cntcnv = 0;
        while (my $line = <I>){
            chomp $line;
            next if ($line =~m/^genesymbol/);
	    next if ($line =~m/^\#N\/A/);
	    next if ($line =~m/^Field4/);

	    my @line = split (/\t/, $line);
	    my $ctg = $line[2];
	    my $start = $line[3];
#	    print STDERR "$start\n";
	    my $end = $line[4];
	    my $pval = $line[6];
	    my $gainloss = $line[7];

	    my $int = $start . "-" . $end;
	    if ($pval eq "NA"){
		next;
	    }
	    elsif ($pval eq ""){
		next;
	    }
	    elsif ($pval > 0.05){
		next;
	    }
	    else{ 
		$cnvs->{$ctg}->{$int}->{$cfile} = $gainloss;
		$cnvclass->{$ctg}->{$int} = "CNV";
		$cntcnv++;
	    }
	}
	close (I);
	print STDERR "$cntcnv\n";
    }
}

# get the unique and exhaustive set of queries
my $qss = {};
foreach my $qs (keys %$qsnvs){
    $qss->{$qs} = 1;
}
foreach my $qc (keys %$qcnvs){
    $qss->{$qc} = 1;
}

# post process datastruc to fill in ref alleles                                                              
# for queries with snp gaps. default to hg ref.                                                               
print STDERR "Post-processing SNPs\n";
if ($snpsdir){
    foreach my $ctg (keys %$snps){
	foreach my $pos (keys %{$snps->{$ctg}}){
	    foreach my $query (keys %$qss){
		if (exists($snps->{$ctg}->{$pos}->{$query})){
		    next;
		}
		else {
#                   if (exists($snps->{$ctg}->{$pos}->{$bloodname})){                                        
#                       $snps->{$ctg}->{$pos}->{$query} = $snps->{$ctg}->{$pos}->{$bloodname};                
#                   }                                                                                        
#                   else{                                                                                     
		    $snps->{$ctg}->{$pos}->{$query} = $refbases->{$ctg}->{$pos};
#                   }                                                                                        
		}
	    }
	}
    }
}

# post process indels if indels are being characterized                                                       
if ($indelsyes){
    print STDERR "Post-processing INDELS\n";
    $pindels = screens ("INDELS", $indels, $class, $qss);
}

if ($cnvsdir){
    print STDERR "Post-processing CNVS\n";
    $pcnvs = screens ("CNVS", $cnvs, $cnvclass, $qss);
}


# collect all the data
my $countermerge = 0;
my @data;
if ($snpsdir){
    push (@data, $snps);
    $countermerge++;
}
if ($indelsyes){
    push (@data, $pindels);
    $countermerge++;
}
if ($cnvsdir){
    push (@data, $pcnvs);
    $countermerge++;
}

# bail if nothing defined
if ($countermerge == 0){
    print STDERR "Select something!\n";
    die;
}

# print the merged matrix
print_matrix (\@data, $class, $cnvclass, $outdir, $qss);



#####SUBS#####

sub screens {
    my $kind = shift;
    my $indels = shift;
    my $class = shift;
    my $qs = shift;
    
    my @queries = keys %$qs;
    
    my $pindels = {};
    foreach my $ctg (keys %$indels){
        foreach my $pos (keys %{$indels->{$ctg}}){
            print STDERR "$ctg\t$pos\n";

            # screen out overlapping elements (most likely imperfect deletions)                              
	    # applicable for CNVs as well
            my $signal = 0;
            my ($posstart, $posend) = split (/\-/, $pos);
            foreach my $int (keys %{$class->{$ctg}}){

                # skip the snps in the interrogation                                                        
                next unless ($int =~m/\-/);

                my ($start, $end) = split (/\-/, $int);

                # skip self                                                                                
                next if (($posstart == $start) and ($posend == $end));

                if (($posstart >= $start) and ($posstart <= $end)){
                    print STDERR "$kind overlap: $class->{$ctg}->{$pos}\t$ctg\t$pos\t$int\n";
                    $signal++;
                    last;
                }
                if (($posend >= $start) and ($posend <= $end)){
                    print STDERR "$kind overlap: $class->{$ctg}->{$pos}\t$ctg\t$pos\t$int\n";
                    $signal++;
                    last;
                }
            }
            next if ($signal > 0);

            # screen out elements that aren't identical (most likely imperfect insertions)                  
	    # not applicable for CNVs where you want to capture a three char state
            my $signal2 = 0;
            my @indels;
            foreach my $query (@queries){
                if (exists($indels->{$ctg}->{$pos}->{$query})){
                    push (@indels, $indels->{$ctg}->{$pos}->{$query});
                }
            }
	    unless (keys %{{ map {$_, 1} @indels }} == 1) {
                $signal2++;
            }
            if (($signal2 > 0) and ($kind eq "INDELS")){
                my $indelprint = join "\t", @indels;
                print STDERR "$kind inconsistent: $class->{$ctg}->{$pos}\t$ctg\t$pos\t$indelprint\n";
                next
            }

            # code the indel if it passes the screens                                                         
	    if ($kind eq "INDELS"){
		foreach my $query (@queries){
		    if (! exists($indels->{$ctg}->{$pos}->{$query})){
			$pindels->{$ctg}->{$pos}->{$query} = 0;
		    }
		    else {
			$pindels->{$ctg}->{$pos}->{$query} = 1;
		    }
		}
	    }

	    # code the cnv if it passes the screens
	    else {
		foreach my $query (@queries){
                    if (! exists($indels->{$ctg}->{$pos}->{$query})){
                        $pindels->{$ctg}->{$pos}->{$query} = 0;
                    }
		    elsif ($indels->{$ctg}->{$pos}->{$query} eq "loss"){
			$pindels->{$ctg}->{$pos}->{$query} = 1;
		    }
		    elsif ($indels->{$ctg}->{$pos}->{$query} eq "gain"){
                        $pindels->{$ctg}->{$pos}->{$query} = 2;
                    }
                    else {
			print STDERR "CNV error\n";
			die;
                    }
                }
            }

            # by definition, the reference does not have either the insertion                                  
            # or deletion or cnv, so assign it a '0'.                                                        
#            $pindels->{$ctg}->{$pos}->{$ref} = 0;
        }
    }
    return ($pindels);
}




sub print_matrix {
    my $data = shift;
    my $class = shift;
    my $cnvclass = shift;
    my $outdir = shift;
    my $qs = shift;

    my @queries = keys %$qs;

    # get some matrix data
    my $taxa = @queries;

    my $alnlentot = 0;
    foreach my $type (@$data){
	foreach my $c (keys %$type){
	    foreach my $s (keys %{$type->{$c}}){
		$alnlentot++;
	    }
	}
    }
    
    # sort print the matrix and charpars    
    print STDERR "Printing matrix\n";                               
    open (CAT, ">$outdir/matrix");
    open (PRT, ">$outdir/partitions");
    print CAT "#NEXUS\n";
    print CAT "BEGIN DATA;\n";
    print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
    print CAT "FORMAT INTERLEAVE SYMBOLS=\"AGCTN012\" MISSING=? GAP=-;\n";
    print CAT "MATRIX\n";
    print PRT "BEGIN SETS;\n";

    my $start = 1;
    my $end;
    foreach my $type (@$data){
	foreach my $contig (keys %$type){
#	    foreach my $count (sort {$a <=> $b} keys %{$snps->{$contig}}){
	    foreach my $count (keys %{$type->{$contig}}){
		$end = $start;

		my $c;
		if ($class->{$contig}->{$count}){
		    $c = $class->{$contig}->{$count};
		}
		elsif ($cnvclass->{$contig}->{$count}){
		    $c = $cnvclass->{$contig}->{$count};
		}
		else {
		    $c = "UNDEF";
		}

		my $countprint = $count;
		$countprint =~s/\-/\_/g;
		
		print CAT "[Partition $contig $c $countprint]\n";
#		my $spcounter = 0;
		foreach my $sp (sort keys %{ $type->{$contig}->{$count} }){
#		    $spcounter++;
#		    my @sp = split (/\_/, $sp);
#		    my $spcount = @sp;
#		    if ($spcount > 1){
#			print CAT "$sp[0]$sp[1]\t$type->{$contig}->{$count}->{$sp}\n";
#		    }
#		    else {
#			print CAT "$sp\t$type->{$contig}->{$count}->{$sp}\n";
#		    }
		    print CAT "$sp\t$type->{$contig}->{$count}->{$sp}\n";
		}
#		print STDERR "$contig $c $countprint $spcounter\n";

		my $charset = $c . "_" . $contig . "_" . $countprint;
		print PRT "CHARSET $charset=$start-$end;\n";
	    
		print CAT "\n";
		$start = $end + 1;
	    }
	}
    }
    print CAT ";\n";
    print CAT "END;\n\n";
    print PRT "END;\n";
    close (CAT);
    close (PRT);
    
    `cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
}

