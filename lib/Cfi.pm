package Cfi;

=head1 NAME

Cfi.pm - contains methods for the Cfi / RPAA pipeline

=head1 SYNOPSIS



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

#####TODO: Document
#####     

use strict;
use warnings;
use Bio::TreeIO;
use Bio::AlignIO;
use Array::Compare;

our $VERSION = "0.01";

sub new {
    my ($class, %args) = @_;
    
    # create the object with room for
    # a tree, and aln, and the config
    my $self = bless {
	'configuration' => undef,
	'alignment'     => undef,
    }, $class;
 
    # if any args, store those
    foreach my $arg (keys %args){
	$self->{$arg} = $args{$arg};
    }
    
    return $self;
}
 
sub load_config{
    my $self = shift;
    my $file = shift;

    # get and store key/val pairs in config
    my $config = {};
    open (F, "$file") or
	die "Unable to open config file\n";
    while (my $line = <F>){
        chomp $line;
	my ($key, $value) = 
	    ($line =~ m/^\s*(\w+)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    $self->{'configuration'} = $config;
}

sub load_aln{
    my $self      = shift;
    my $alignfile = shift;

    open (NEX, "$alignfile");
    my $charset = {};
    my $nexus   = {};
    while (my $line = <NEX>){
	chomp $line;

        # take only first instances of all of these things                                     
        # header only                                                                 
        ($nexus->{'nchar'}    = $1) if (($line =~m/nchar\s*=\s*(\d+)/i) and (!$nexus->{'nchar'}));
        ($nexus->{'ntax'}     = $1) if (($line =~m/ntax\s*=\s*(\d+)/i) and (!$nexus->{'ntax'}));
        ($nexus->{'datatype'} = $1) if (($line =~m/datatype\s*=\s*(\w+)/i) and (!$nexus->{'datatype'}));
        ($nexus->{'missing'}  = $1) if (($line =~m/missing\s*=\s*(.{1})/i) and (!$nexus->{'missing'}));
        ($nexus->{'gap'}      = $1) if (($line =~m/gap\s*=\s*(.{1})/i) and (!$nexus->{'gap'}));

	if ($line =~m/outgroup/i){
	    $line =~s/outgroup//ig;
	    $line =~s/\s+//g;
	    $line =~s/\;//g;
	    
	    # any instances of more than one outgroup???? <====FIX
	    $nexus->{'outgroup'} = $line; 
	}
	
	if ($line =~m/^\s*charset/i){
	    $line =~s/charset//ig;
	    $line =~s/\s+//g;
	    $line =~s/\;//g;
	    
            my ($partition, $coords) =
                split (/\=/, $line);
	    $charset->{$partition} = $coords;
	}
	$nexus->{'charset'} = $charset;
    }
    close (NEX);
    
    $self->{'alignment'} = $nexus;
}

sub parse_tree {
    my $self     = shift;
    my $treefile = shift;
    my $root     = shift;
#    my $treetype = shift;

    my $input = Bio::TreeIO->new(-file      => $treefile,
				 -format    => "nexus");
#                                 -format    => $treetype);

    # now remove the taxa and analyze the culled trees                            
#    my $counter = 0;
#    my $topos   = {};

#    while (my $tree = $input->next_tree()){
    my $tree = $input->next_tree();
    
    # re-root for nodal comparisons, if required
    if ($root){
        my @nodes = $tree->find_node(-id => $root);
        my $node = $nodes[0];
        $tree->reroot($node);
    }
#    $counter++;

    # get all nodes and internal nodes                                
    my @allnodes      = $tree->get_nodes;
    my @internalnodes = grep {!$_->is_Leaf} @allnodes;
    
    my $lineages = {};
    my $internal_counter = 0;
    foreach my $node (@internalnodes){
	$internal_counter++;
	my @descendents = $node->get_all_Descendents;
	my @children    = grep{$_->is_Leaf} @descendents;
	
	my @array;
   	foreach my $leaf (@children){
	    my $id = $leaf->id;
	    push (@array, $id);
	}
	
	# this checks to see that all nodes have a uniq                       
	# list of children --> compensates for a pruning bug                  
	# in TreeIO where an extra set of parens is left                      
	# on outermost nest;                                                  
	my @sortedarray = sort @array;
	my $arrayjoin = join(",", @sortedarray);
	$lineages->{$arrayjoin} = [@sortedarray];
	
    }
    
#    $topos->{$counter} = $lineages;
#}
    
    return ($lineages);
}


sub compare_all_nodes {
    my $self   = shift;
    my $TEtree  = shift;
    my $CFItree = shift;

    my $topo = {};
    my $treesame = 1; #assume the tree is the same topo; if find one node off, report that it is not
    foreach my $nodeTE (sort keys %$TEtree){
        my @sortedleavesTE = sort @{ $TEtree->{$nodeTE} };
        my $sortedleavesTE = join (";", @sortedleavesTE);
	
	my $match = 0;
	foreach my $nodeCFI (sort keys %$CFItree){
	    my @sortedleavesCFI = sort @{ $CFItree->{$nodeCFI} };
	    my $sortedleavesCFI = join (";", @sortedleavesCFI);
	    
	    my $comp = Array::Compare->new;
	    if ($comp->compare(\@sortedleavesTE, \@sortedleavesCFI)) {
		$match = 1;
		last;
	    }
	    else {
		next;
	    }
	}
	
	if ($match == 1){
	    $topo->{$nodeTE} = 1;
	}
	else {
	    $topo->{$nodeTE} = 0;
	    $treesame = 0;
	}
    }
    return ($topo, $treesame);
}


sub integrate_trap{
    my $self = shift;
    my $x = shift;
    my $y = shift;

    my $points = @$x;
    my $integration = 0;
    for (my $i = 0; $i < $points - 1; $i++){
        my $height = @$x[$i+1] - @$x[$i];
        my $area   = 0.5 * $height * (@$y[$i+1] + @$y[$i]);
        $integration += $area;
#	print STDERR "$i\t$height\t$area\t$integration\n";
    }

    return ($integration);
}

sub parse_loess{
    my $self = shift;
    my $infile = shift;

    my @x;
    my @y;
    open (FILE, "$infile");
    while (my $line = <FILE>){
	chomp $line;
	my ($x, $y) = split (/\t/, $line);
	push (@x, $x);
	push (@y, $y);
    }
    close (FILE);

    return (\@x, \@y);
}


sub get_charsets{
    my $self = shift;
    return $self->{'alignment'}->{'charset'};
}

sub get_ntax{
    my $self = shift;
    return $self->{'alignment'}->{'ntax'};
}

sub get_lineages{
    my $self = shift;
    return $self->{'tree'};
}

sub get_root{
    my $self = shift;
    return $self->{'configuration'}->{'ROOT'};
}

sub get_mode{
    my $self = shift;
    return $self->{'configuration'}->{'MODE'};
}

sub run_paup{
    my $self   = shift;
    my $file   = shift;    
    `paup -n $file`;
}

sub run_raxml{
    my $self = shift;
    my $file = shift;
    my $root = shift;
    my $outdir = shift;
    
    my @filename = split (/\//, $file);
    my $filename = pop @filename;

    my $cmd = "raxmlHPC-SSE3 ";
    $cmd .= "$self->{'configuration'}->{'TREECOMMAND'} ";
    $cmd .= "-o $root ";
    $cmd .= "-s $file ";
    $cmd .= "-n $filename";
    
    `$cmd`;
    `mv RAxML_info.$filename $outdir/logs`;
    `mv RAxML_log.$filename $outdir/logs`;
    `mv RAxML_bestTree.$filename $outdir/trees`;
    `rm RAxML_parsimonyTree.$filename`;
    `rm RAxML_result.$filename`;
}

sub generate_consensus{
    my $self = shift;
    my $outdir = shift;
    my $matrixfile = shift;
    my $treesref = shift;
    my $name = shift;
    
    # head                                                                                               
    open (CFI, ">$outdir/consensus/$name.nex");
    print CFI "#nexus;\n";
    print CFI "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
    print CFI "set increase=auto;\n";
    print CFI "log start file=$outdir/consensus/$name.out replace=yes;\n";
    print CFI "execute $matrixfile;\n";
    
    my $firsttree = @$treesref[0];
    print CFI "gettrees file=$outdir/trees/$firsttree mode=3;\n";

    foreach my $tree (@$treesref){
	next if $tree eq $firsttree;
	print CFI "gettrees file=$outdir/trees/$tree mode=7;\n";
    }
    
    print CFI "outgroup $self->{'configuration'}->{'ROOT'};\n";
    print CFI "$self->{'configuration'}->{'CONTREECOMMAND'} treefile=$outdir/consensus/$name.tre;\n";
    close (CFI);
}

sub get_rscript{
    my $self = shift;
    return ($self->{'configuration'}->{'RSCRIPT'});
}
 
sub get_span{
    my $self = shift;
    return ($self->{'configuration'}->{'SPAN'});
}

sub get_sections{
    my $self = shift;
    return ($self->{'configuration'}->{'SECTIONS'});
}

sub round {
    my($number) = shift;
    return int($number + .5);
}


sub generate_cmds{
    my $self = shift;
    my $mode = shift;
    my $charset    = shift;
    my $maxnum     = shift;
    my $step       = shift;
    my $outdir     = shift;
    my $matrixfile = shift;

    my $sampleMaster = {};

    # generate the steps
    my %dothese;
    my $counter = 0;
    while (1){
	$counter = $counter + $step;
	last if ($counter > $maxnum);
	$dothese{$counter} = 1;
    }
    
    # always do the first and second cases
    ($dothese{'1'} = 1) unless (exists ($dothese{'1'}));
    ($dothese{$maxnum} = 1) unless (exists ($dothese{$maxnum}));
    my @dothese = sort {$a <=> $b} keys %dothese;
    
    for (my $sample = 1; $sample <= $self->{'configuration'}->{'SAMPLES'}; $sample++){
        my %charset_cache = %{$charset};
        my @includes;

	foreach my $number (@dothese){

	    # include all partitions even if stepping around!
	    while (scalar (@includes) < $number){
		my @partitions = keys %charset_cache;
		my $key        = int rand(@partitions);
		push (@includes, $partitions[$key]);
		delete $charset_cache{$partitions[$key]};
	    }

	    my $includes_string = join(" ", @includes);
	    $sampleMaster->{$sample}->{$number} = $includes_string;
	    
	    # skip the rest if we're in ML mode
	    if ($mode eq "ML"){
		open (CFI, ">$outdir/cmds/$sample.$number.cmds");
		print CFI "$includes_string\n";
		close (CFI);
	    }
	    elsif ($mode eq "PARSIMONY"){
	    
		# head                                                                                   
		open (CFI, ">$outdir/cmds/$sample.$number.cmds");
		print CFI "#nexus;\n";
		print CFI "set warntree=no warntsave=no warnreset=no notifybeep = no;\n";
		
		# maxtrees and logs                                                                
		if ($self->{'configuration'}->{'MAXTREES'}){
		    print CFI "set maxtrees=$self->{'configuration'}->{'MAXTREES'} increase=no;\n";
		}
		else {
		    print CFI "set increase=no;\n";
		}
		
		print CFI "log start file=$outdir/logs/$sample.$number.out replace=yes;\n";
		
		# search                                                                                 
		print CFI "execute $matrixfile;\n";
#		print CFI "pset gapmode=newstate;\n"; # <--CUSTOM CODE!!
		print CFI "exclude all;\n";
		
		my $list = join " ", @includes;
		print CFI "include $list;\n";
		print CFI "$self->{'configuration'}->{'TREECOMMAND'};\n";
		
		# CONTREE here to create consensus multiple 
		# most parsimonious trees if exists
		print CFI "contree all/strict=yes majrule=no showtree=no treefile=$outdir/trees/$sample.$number.tre;\n";
		print CFI "cleartrees;\n";
		print CFI "log stop;\n";
		print CFI "quit /warnTsave=no;\n";
		close (CFI);
	    }
	    else {
		print STDERR "Unknown run mode\n";
		die;
	    }
	}
    }	
 
   return $sampleMaster;
}

sub sub_matrix{
    my $self = shift;
    my $sample = shift;
    my $number = shift;
    my $parties = shift;
    my $matrixfile = shift;
    my $charsets = shift;
    my $outdir = shift;

    # output sectioned matrix
    my $alnin = Bio::AlignIO->new(-file => "$matrixfile", -format => "nexus");
    my $alnobj = $alnin->next_aln();
    
    my $alndata = {};
    my $alnlens = {};
    my $alnlentot = {};
    my %taxa;
    
    my @parties = split (/ /, $parties);
    foreach my $seq ($alnobj->each_seq){
	my $id        = $seq->display_id;
	$taxa{$id} = 1;
	
	foreach my $party (@parties){
	    my $coords = $charsets->{$party};
	    my ($start, $end) = split (/\-/, $coords);
	    my $partition = $seq->subseq($start, $end);
	    my $alnlen = $end - $start + 1;
	    $alnlentot->{$id} += $alnlen; #all will be the same                                  
	    $alnlens->{$party} = $alnlen;
	    
	    $alndata->{$party}->{$id} = $partition;
	}
    }
    
    # sort print the concatenation and charpars                                                 
    my @taxa = keys %taxa;
    my $taxa = @taxa;
    open (CAT, ">$outdir/$sample.$number.matrix.tmp");
    open (PRT, ">$outdir/$sample.$number.partitions.tmp");
    print CAT "#NEXUS\n";
    print CAT "BEGIN DATA;\n";
    print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot->{$taxa[0]};\n";
    print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
    print CAT "MATRIX\n";
    print PRT "BEGIN SETS;\n";
    
    my $start = 1;
    my $end;
    foreach my $count (sort keys %$alndata){
	$end = $start - 1 + $alnlens->{$count};
	print CAT "[Partition $count length $alnlens->{$count} chars $start-$end]\n";
	print PRT "CHARSET $count=$start-$end;\n";
	foreach my $sp (sort keys %{ $alndata->{$count} }){
	    print CAT "$sp\t$alndata->{$count}->{$sp}\n";
	}
	print CAT "\n";
	$start = $end + 1;
    }
    
    print CAT ";\n";
    print CAT "END;\n\n";
    print PRT "END;\n";
    close (CAT);
    close (PRT);
    
    `cat $outdir/$sample.$number.matrix.tmp $outdir/$sample.$number.partitions.tmp > $outdir/$sample.$number.nex`;
    `rm $outdir/$sample.$number.matrix.tmp`;
    `rm $outdir/$sample.$number.partitions.tmp`;
}


sub aln_converter{
    my $self = shift;
    my $infile = shift;
    my $informat = shift;
    my $outfile = shift;
    my $outformat = shift;
    
    my $alnin = Bio::AlignIO->new(-file   => "$infile",
				  -format => "$informat");
    my $alnout = Bio::AlignIO->new(-file   => ">$outfile",
				   -format => "$outformat");
    
    while (my $aln = $alnin->next_aln){
	$alnout->write_aln($aln);
    }
}

sub tree_converter{
    my $self = shift;
    my $infile = shift;
    my $informat = shift;
    my $outfile = shift;
    my $outformat = shift;
    
    my $treein = new Bio::TreeIO(-file => "$infile",
				 -format => "$informat");
    my $treeout = new Bio::TreeIO(-file => ">$outfile",
				  -format => "$outformat");
    while (my $tree = $treein->next_tree){
	$treeout->write_tree($tree);
    }
    
    # what follows is an intense hack necessary because the bioperl 
    # modules wrap the nexus tree in an extra set of parens
    open (TREE, "$outfile");
    open (TREEO, ">$outfile.tmp");
    while (my $line = <TREE>){
	chomp $line;
	if ($line =~m/tree\sBioperl\_/){
	    $line =~s/\[\]//g;
	    my @line = split (/\s+/, $line);
	    $line[5] =~s/^.//;
	    $line[5] =~s/.{2}$//;
	    $line[5] .= ";";
	    print TREEO "tree Bioperl_1 = [&U] $line[5]\n";
	}
	else {
	    print TREEO "$line\n";
	    next;
	}
    }

    `rm $outfile`;
    `mv $outfile.tmp $outfile`;

    close (TREE);
    close (TREEO);
}

1;
