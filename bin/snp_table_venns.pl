#!/usr/bin/perl

# -i is the tabe from snp_matrixBuilder.pl
# then input multiple arrays of strains in @ARGV

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $in    = $opts{'i'};

# store all the strain lists and tie them together
my $groups = {};
foreach my $file (@ARGV){
    open (F, "$file");
    my $file_name;
    if ($file =~/\//g){
        $file =~m/.*\/(.*)$/;
        $file_name = $1;
    }
    else {
        $file_name = $file;
    }

    while (my $line = <F>){
	chomp $line;
#	$groups->{$file_name} = $line;
	$groups->{$line} = $file_name;
    }
}

# go through the snp table line by line
my $venn = {};
my $venncount = {};
my $strains = {};
open (I, "$in");
while (my $line = <I>){
    chomp $line;

    # index the strains
    if ($line =~m/Contig\tSNP/){
	my @data2 = split (/\t/, $line);
	shift @data2;
	shift @data2;
	my $counter = 0;
	foreach my $strain (@data2){
	    $counter++;
	    $strains->{$counter} = $strain;
	}
	next;
    }

    # sift the columns and mine them for position and reference
    # assumes that the data has annotation information!!!
    my @data0 = split (/\t/, $line);
    pop (@data0); pop (@data0); pop (@data0); pop (@data0); pop (@data0); pop (@data0); pop (@data0);
    shift (@data0);
    my $pos = shift (@data0);
    my $ref = pop (@data0);

    # assume that ambiguous characters are not snps and that they should
    # have the identify of the reference
    my @data;
    foreach my $piece (@data0){
	if ($piece eq "?"){
	    push (@data, $ref);
	}
	else {
	    push (@data, $piece);
	}
    }
    
    # check to see if all the strains share the same snp
    if (keys %{{ map {$_, 1} @data }} == 1){
	$venn->{'ALL'}->{$pos} = 1;
	$venncount->{'ALL'}++;
	next;
    }

    # otherwise cycle through the strains and sort strain groups
    # by unique alleles they contain
    else{
	my $int = {};
	my $cnt = 0;

	# store alleles by group
	foreach my $s (@data){
	    $cnt++;
	    my $str = $strains->{$cnt};
	    $int->{$groups->{$str}}->{$s}++;
	}
	
	# look at PE specific alleles
	my $pe = $int->{'list.pe'};
#	my @peall;
#	foreach my $p (sort keys %$pe){
#	    my $pcnt = $pe->{$p};
#	    my $pstr = "$p($pcnt)";
#	    push (@pstr, @peall);
	my @pestring = sort keys %$pe;
	my $pestring = join ",", @pestring;

	my $se2 = $int->{'list.se2'};
        my @se2string = sort keys %$se2;
        my $se2string = join ",", @se2string;
	
	my $se = $int->{'list.se_minusbad'};
        my @sestring = sort keys %$se;
        my $sestring = join ",", @sestring;
	
	if (($pestring ne $se2string) and ($pestring ne $sestring)){
	    print "NC_000962\t$pos\t$pos\n";
	}
    }
}
=head	
	# invert and store groups by unique allele string
	my $invertint = {};
	foreach my $grp (sort keys %$int){
	    my @residues = sort keys %{$int->{$grp}};
            my $residues = join ",", @residues;
	    $invertint->{$residues}->{$grp} = 1;
	}
	
	my @reportstring;
	foreach my $allele (sort keys %$invertint){
	    my @grpstring;
	    foreach my $grp (sort keys %{$invertint->{$allele}}){
		push (@grpstring, $grp)
	    }
	    my $grpstring = join ",", @grpstring;
	    push (@reportstring, $grpstring);
	}
	my $reportstring = join "\t", @reportstring;
	print "$pos\t$reportstring\n";
    }
}	    
	print yellow;
#	my @key;
#	foreach my $grp1 (sort keys %$int){
#	    my @grp1;
#	    my @residues1 = sort keys %{$int->{$grp1}};
#	    my $residues1 = join ",", @residues1;
#	    foreach my $grp2 (sort keys %$int){
#		next if ($grp2 eq $grp1);
#		my @residues2 = sort keys %{$int->{$grp2}};
#		my $residues2 = join ",", @residues2;

#		if ($residues1 eq $residues2){
#		    push (@grp1, $grp2);
		    

#    }
#}








    my $rate;
#    if ($data[2] eq $data[3] eq $data[4] eq $data[5]){
#    if ($data[2] eq $data[3] eq $data[4]){
    if (keys %{{ map {$_, 1} @data2 }} == 1) {
	$rate = "N";
    }
    else {
	if (($data[2] eq "?") or ($data[3] eq "?") or ($data[4] eq "?") or ($data[5] eq "?")){
	    $rate = "+";
	}
	else {
	    $rate = "*";
	}
    }
    
    print "$line\t$rate\n";
    
}
close (I);
=cut
