package TreeSupports;

=head1 NAME

TreeSupports.pm - contains methods for caculating BS, PBS, PHBS, and HS

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2012 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####TODO:                                                                        
#####      Fix for multiple pbs trees
#####      Rooting issue when pruning -- lose the true root?
#####          current fix is a hack --> uniqifying the sub nodes
#####      Will lineages always be consistent across the three phases?
#####      Also -- deal with multiple outgroups??                                 
#####      Multiple gene partitions                                               
#####      Should always be analyzing nodes = tax - 3 (we have one extra)         
#####      Add function for consensus trees            
#####      Documentation and Formatting

use strict;
use warnings;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::SeqIO;
use Array::Compare;

our $VERSION = "0.01";

sub new {
    my ($class, %args) = @_;
    
    # create the object with room for
    # a tree, and aln, and the config
    my $self = bless {
	'tree'          => undef,
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

    # check for PAUP in user's path if not in config                                                 
    if (! exists($config->{'PAUP'})){
	my $paup = `which paup`;
	chomp $paup;

        if (defined ($paup)){
            $config->{'PAUP'} = $paup;
        }
        else {
            die "No access to PAUP!\n";
	}
    }

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
	
        if ($line =~m/\s*charset/i){
            my ($partition, $coords) =
                split (/\=/, $line);

            $partition =~s/charset\s*//ig;
            $partition =~s/\s+//g;

            $coords =~s/\;//g;
            $coords =~s/^\s//;

            $charset->{$partition} = $coords;
	}

	$nexus->{'charset'} = $charset;
    }
    close (NEX);
    
    $self->{'alignment'} = $nexus;
}

sub load_tree{
    my $self     = shift;
    my $treefile = shift;
    my $root     = shift;
    my $type     = shift;

    # generate the tree object                                                    
    my $in_tree  = Bio::TreeIO->new(-file   => $treefile,
                                    -format => $type);
    my $tree     = $in_tree->next_tree;

    # re-root if required
    if ($root){
	my @nodes = $tree->find_node(-id => $root);
	my $node = $nodes[0];
	$tree->reroot($node);
    }

    # get all nodes and internal nodes                                            
    my @allnodes      = $tree->get_nodes;
    my @internalnodes = grep {!$_->is_Leaf} @allnodes;

    # cycle through internal nodes and get each internal                          
    # node's immediate descendents. Get all leaves from                           
    # each immediate descendent and store.                                        
    my $lineages = {};
    my $internal_counter = 0;
    foreach my $node (@internalnodes){
        $internal_counter++;
	my @lineages = $node->each_Descendent;

        my $lineage_counter = 0;
        foreach my $lineage (@lineages){
            $lineage_counter++;

            my @children;
            if ($lineage->is_Leaf){
                push (@children, $lineage);
            }
            else {
                my @descendents = $lineage->get_all_Descendents;
                @children = grep{$_->is_Leaf} @descendents;
            }

            foreach my $leaf (@children){
                my $id = $leaf->id;
		push @{$lineages->{$internal_counter}->{$lineage_counter}}, $id;
            }
        }
    }
    $self->{'tree'} = $lineages;

}

sub prune_tree {
    my $self     = shift;
    my $taxa     = shift;
    my $treefile = shift;

    my $input = Bio::TreeIO->new(-file      => $treefile,
				 -format    => "nexus");

    # now remove the taxa and analyze the culled trees
    my $counter = 0;
    my $topos   = {};
    while (my $tree = $input->next_tree()){
	$counter++;
	foreach my $taxon (keys %$taxa) {
	    $tree->remove_Node($taxon);
	}

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
#		push @{$lineages->{$internal_counter}}, $id;
	    }
	
	    # this checks to see that all nodes have a uniq
	    # list of children --> compensates for a pruning bug
	    # in TreeIO where an extra set of parens is left
	    # on outermost nest;
	    my @sortedarray = sort @array;
	    my $arrayjoin = join(":", @sortedarray);
	    $lineages->{$arrayjoin} = [@sortedarray];
#	    my $arrayjoin = join(":", @array);
#            $lineages->{$arrayjoin} = [@array];

	}
	
	$topos->{$counter} = $lineages;
    }
    
    return ($topos);
}

sub compare_trees {
    my $self   = shift;
    my $tree1  = shift;
    my $tree2  = shift;
    
    my @leafstring;
    foreach my $node (sort keys %$tree1){
	my @sortedleaves = sort @{ $tree1->{$node} };
	my $sortedleaves = join (";", @sortedleaves);
	push (@leafstring, $sortedleaves);
    }

    my @leafstringTE;
    foreach my $nodeTE (sort keys %$tree2){
	my @sortedleavesTE = sort @{ $tree2->{$nodeTE} };
	my $sortedleavesTE = join (";", @sortedleavesTE);
	push (@leafstringTE, $sortedleavesTE);
    }

    my @sortedleafstring   = sort @leafstring;
    my @sortedleafstringTE = sort @leafstringTE;
    my $comp = Array::Compare->new;

    if ($comp->compare(\@sortedleafstring, \@sortedleafstringTE)) {
	return ("same", \@sortedleafstring, \@sortedleafstringTE);
    }
    else{
	return ("diff", \@sortedleafstring, \@sortedleafstringTE);
    }
}
    

sub transform_aln{
    my $self      = shift;
    my $alignfile = shift;
    my $informat  = shift;
    my $outformat = shift;
    my $outdir    = shift;
    my $charset   = shift;
    my $coords    = shift;

    my $alnin = Bio::AlignIO->new(-file   => "$alignfile",
				  -format => "$informat");
        
    # if the coords are a range, do this:                                                                    
#    if ($coords =~m/\-/){
    my ($start, $end) = split (/\-/, $coords);
	
    open (FA, ">$outdir/$charset.$outformat");
    while (my $aln = $alnin->next_aln()){
	foreach my $seq ($aln->each_seq){
	    my $id        = $seq->display_id;
	    my $partition = $seq->subseq($start, $end);
	    
	    print FA ">$id\n$partition\n";
	}
    }
    close (FA);
#    }	

    # otherwise they are a list of positions                                                                 
#    else {
#	my @coords = split (/\s/, $coords);
#
#	open (FA, ">$outdir/$charset.$outformat");
#	foreach my $coord (@coords){
#	    while (my $aln = $alnin->next_aln()){
#		foreach my $seq ($aln->each_seq){
#		    my $id        = $seq->display_id;
#		    my $partition = $seq->subseq($coord, $coord);
#
#		    print FA ">$id\n$partition\n";
#		}
#	    }
#	}
#	close (FA);
#    }
}

sub find_inf_taxa{
    my $self      = shift;
    my $alignfile = shift;
    my $outformat = shift;
    my $outdir    = shift;
    my $cset      = shift;

    my $missing = $self->{'alignment'}->{'missing'};
   
    my %full_taxa;
    my %all_taxa;
    my %missing_taxa;
    my $seqin = Bio::SeqIO->new(-file   =>"$outdir/$cset.$outformat",
				-format =>"$outformat");
    while (my $seq = $seqin->next_seq()){
        my $id             = $seq->display_id;
        my ($acc, $coords) = split (/\//, $id);
        my $length         = $seq->length;
	my $sequence       = $seq->seq;

        my $matchstring = '\\' . $missing . "{$length}"; # a bit hokey!           
        ($full_taxa{$acc} = 1) unless ($sequence =~m/$matchstring/);
	($missing_taxa{$acc} = 1) if ($sequence =~m/$matchstring/);
        $all_taxa{$acc} = 1;

    }
    return (\%full_taxa, \%all_taxa, \%missing_taxa);
}

sub check_topo_build_bs{
    my $self        = shift;
    my $alignfile   = shift;
    my $treefile    = shift;
    my $full_taxa   = shift;
    my $all_taxa    = shift;
    my $missing_taxa= shift;
    my $outdir_cmds = shift;
    my $outdir_logs = shift;
    my $outdir_prun = shift;
    my $indir_pbs   = shift;
    my $cset        = shift;
    my $coords      = shift;
    my $test        = shift;

    my $lineages = $self->{'tree'};
    my $config   = $self->{'configuration'};

    open (BS, ">$outdir_cmds/$cset.bs_commands.nex");

    # print out the header                                                        
    print BS "#NEXUS\n";
    print BS "set increase=auto;\n";
    print BS "set WarnRoot=no;\n";

    # execute the nexus aln file                                                  
    print BS "execute $alignfile;\n";

    # cycle through nodes, lineages and leaves to see what has data               
    # and cache in lineage specific vars.
    # check to see if the topology of the                                                                    
    # PBS without tree is the same as the topology of                                          
    # the SA tree with respect to informative taxa.                                             
    # Write the information.          

    my @savednodes;
    foreach my $node (sort keys %$lineages){

	# check to see if the topologies of the pruned trees are the same
	# coded to compare the without tree only
	open (CMP, ">$outdir_prun/$cset.$node.prun");
	if (-e "$indir_pbs/$node.without.tre"){
	    my $topospbs = $self->prune_tree ($missing_taxa, "$indir_pbs/$node.without.tre");
	    my $toposTE  = $self->prune_tree ($missing_taxa, $treefile);

	    foreach my $tree2 (sort keys %$topospbs){
		foreach my $tree1 (sort keys %$toposTE){
		    
		    my ($compare, $tpbs, $tTE) = 
			$self->compare_trees ($topospbs->{$tree2}, $toposTE->{$tree1});
		    print CMP "$compare\t";
		    my $pbs_string = join (":", @$tpbs);
		    my $te_string  = join (":", @$tTE);
		    print CMP "$pbs_string\t$te_string\n";
		}
	    }
	}
	else {
	    print CMP "skipped\n";
	}
	close (CMP);
		
	
        # check to see if this node has information                                             
        # in each of its lineages; should also deal with polytomies 
	
	# check to see that there are at least two defined                                             
        # taxa outside of our constraint (one required for                                           
        # BS relevance and 2 required for PAUP to run).                                            
        # this also disallows analysis of the root node 
	
	# check to see if we are enforcing only two existing taxa
	# across all lineages, or one to exist in each lineage
        my %constraintleaves;
        my %lineage_tracker;
	my @lineage_cache;
	my $simple_taxa_counter = 0;

        foreach my $lineage (sort keys %{ $lineages->{$node} }){
            push (@lineage_cache, $lineage);

            foreach my $leaf (@{ $lineages->{$node}->{$lineage} }){

		$constraintleaves{$leaf} = 1;

                if (exists ($full_taxa->{$leaf})){
                    $lineage_tracker{$lineage} = 1;
		    $simple_taxa_counter++;
                }
                else {
		    next;
                }

            }
	}

        my @lineage_tracker = keys %lineage_tracker;
        my $lineage_tracker = @lineage_tracker;
        my $lineage_cache   = @lineage_cache;

        my $noncon_defined = 0;
        my @nonconstraintleaves;
        foreach my $taxa (sort keys %$all_taxa){
            (push (@nonconstraintleaves, $taxa))
                unless (exists ($constraintleaves{$taxa}));
        }

        foreach my $noncon (@nonconstraintleaves){
            ($noncon_defined++) if (exists ($full_taxa->{$noncon}));
        }

	# check to see what we are testing in terms of nodes
        # if conditions met, print out the constraint and store;                  
        # if not, iterate to the next node                                        
	if ($test == 1){
	    if (($lineage_tracker == $lineage_cache) and ($noncon_defined >= 2)){
		print BS "CONSTRAINTS $node=((";
		push (@savednodes, "$node");
		
		print BS join (",", sort keys %constraintleaves);
		print BS "));\n";
	    }
	    else {
		next;
	    }
	}
	else {
	    if (($simple_taxa_counter >= 2) and ($noncon_defined >= 2)){
		print BS "CONSTRAINTS $node=((";
		push (@savednodes, "$node");

                print BS join (",", sort keys %constraintleaves);
                print BS "));\n";
	    }
	    else {
		next;
	    }
	}
    }    
    
    # print out the taxset                                                        
    print BS "Begin sets;\n";
    print BS "TAXSET informative=";
    print BS join (" ", sort keys %$full_taxa);
    print BS ";\n";
    print BS "end;\n";
    
    # for this charset calculate the BS values for each                           
    # informative node after restricting to informative taxa                      
    print BS "[!==>BS for charset $cset]\n";
    print BS "log start file=$outdir_logs/$cset.bs.log replace=yes;\n";
    print BS "exclude all;\n";
    print BS "include $coords;\n";
    print BS "delete all;\n";
    print BS "restore informative;\n";
    
    foreach my $node (@savednodes){
	print BS "[!==>constraint: $node]\n";
	
	for my $state ("yes", "no"){
	    print BS "$config->{TREECOMMAND} ";
	    print BS "enforce converse=$state constraints=$node;\n";
	    print BS "lenfit;\n";
	    print BS "[!==>end_record]\n";
	}
    }
    
    print BS "log stop;\n";
    print BS "quit /warnTsave=no;\n";
    close (BS);
    
    return ($cset);
}

sub build_ndi{
    my $self       = shift;
    my $alignfile  = shift;
    my $node       = shift;
    my $outdir_cmds = shift;
    my $outdir_logs = shift;
    my $outdir_trees = shift;
    my $enforce = shift;

    my $lineages = $self->{'tree'};
    my $config   = $self->{'configuration'};
    my $charsets = $self->{'alignment'}->{'charset'};
    my $nchar    = $self->{'alignment'}->{'nchar'};

    open (S, ">$outdir_cmds/$node.$enforce.ndi_commands.nex");

    # print out the header                                                                       
    print S "#NEXUS\n";
    print S "set increase=auto;\n";
    print S "set WarnRoot=no;\n";

    # execute the nexus aln file                                                                   
    print S "execute $alignfile;\n";

    # print out the constraint for this node                                                                 
    my @lins;
    foreach my $lin (sort keys %{ $lineages->{$node} }){
	push (@lins, @{ $lineages->{$node}->{$lin} });
    }

    print S "CONSTRAINTS $node=((";
    print S join (",", @lins);
    print S "));\n";
    
    # PUAP setup                                                                                         
    print S "[!==>NDI for node $node]\n";
    print S "log start file=$outdir_logs/$node.$enforce.ndi.log replace=yes;\n";

    # all partitions calculation                                                                   
    print S "[!==>partition: all]\n";
    print S "include all;\n";

    # one search on ALL the data!                                                          
    # include and exclude from there                                                          
    print S "$config->{TREECOMMAND} ";

    # constraint direction                                                                                   
    if ($enforce eq "without"){
	print S "enforce converse=yes constraints=$node;\n";
    }
    if ($enforce eq "with"){
	print S "enforce converse=no constraints=$node;\n";
    }
    
    print S "filter best;\n";
    print S "savetrees file=$outdir_trees/$node.$enforce.tre format=nexus replace=yes root=yes;\n";
    print S "[!==>end_record]\n";

    # searches after charset removals, no trees kept
    foreach my $charset (keys %$charsets){

	# if the coords are a range, do this:                                                           
	my $partLen;
	if ($charsets->{$charset} =~m/\-/){
	    my ($start, $end) = split (/\-/, $charsets->{$charset});
	    $partLen       = $end - $start + 1;
	}

	# otherwise they are a list of positions                                                     
	else {
	    my @coords = split (/\s/, $charsets->{$charset});
	    $partLen = @coords;
	}

	# skip if you're going to exclude all chars
	next if ($partLen >= $nchar);
	
        print S "[!==>partition: $charset]\n";
        print S "include all;\n";
        print S "exclude $charsets->{$charset};\n";
	print S "$config->{TREECOMMAND} ";

	if ($enforce eq "without"){
	    print S "enforce converse=yes constraints=$node;\n";
	}
	if ($enforce eq "with"){
	    print S "enforce converse=no constraints=$node;\n";
	}

        print S "[!==>end_record]\n";
    }

    print S "log stop;\n";
    print S "quit /warnTsave=no;\n";
    close (S);

    return ($node);
}


sub build_pbs{
    my $self       = shift;
    my $alignfile  = shift;
    my $node       = shift;
    my $outdir_cmds = shift;
    my $outdir_logs = shift;
    my $outdir_trees = shift;
    my $enforce     = shift;

    my $lineages = $self->{'tree'};
    my $config   = $self->{'configuration'};
    my $charsets = $self->{'alignment'}->{'charset'};

    open (S, ">$outdir_cmds/$node.$enforce.pbs_commands.nex");

    # print out the header                                                        
    print S "#NEXUS\n";
    print S "set increase=auto;\n";
    print S "set WarnRoot=no;\n";

    # execute the nexus aln file                                                  
    print S "execute $alignfile;\n";

    unless ($node eq "unconstrained"){

        # print out the constraint for this node                                  
        my @lins;
        foreach my $lin (sort keys %{ $lineages->{$node} }){
            push (@lins, @{ $lineages->{$node}->{$lin} });
	}

        print S "CONSTRAINTS $node=((";
        print S join (",", @lins);
        print S "));\n";

    }

    # PUAP setup                                                                  
    print S "[!==>PBS for node $node]\n";
    print S "log start file=$outdir_logs/$node.$enforce.pbs.log replace=yes;\n";

    # all partitions calculation                                                  
    print S "[!==>partition: all]\n";
    print S "include all;\n";

    # only one search on ALL the data!                                            
    # include and exclude from there                                              
    if ($node eq "unconstrained"){
        print S "$config->{TREECOMMAND};\n";
    }
    else {
        print S "$config->{TREECOMMAND} ";
	
	# constraint direction
	if ($enforce eq "without"){
	    print S "enforce converse=yes constraints=$node;\n";
	}
	if ($enforce eq "with"){
	    print S "enforce converse=no constraints=$node;\n";
	}

    }
    print S "filter best;\n";
    print S "savetrees file=$outdir_trees/$node.$enforce.tre format=nexus replace=yes root=yes;\n";
    print S "[!==>end_record]\n";
    
    foreach my $charset (keys %$charsets){
        print S "[!==>partition: $charset]\n";
        print S "exclude all;\n";
        print S "include $charsets->{$charset};\n";
        print S "lenfit;\n";
        print S "[!==>end_record]\n";
    }

    print S "log stop;\n";
    print S "quit /warnTsave=no;\n";
    close (S);
    
    return ($node);
}

sub parse_bs{
    my $self   = shift;
    my $file   = shift;

    my $constraint;
    my $parsinf;
    my $numtrees;
    my $conv;
    my $length;
    
    my %bsdata;
    open (BS, "$file");
    while (my $line = <BS>){
        chomp $line;
	
        # mine the data we need                                                   
        ($constraint = $1) if ($line =~m/==>constraint:\s*(\w+)/i);
	($parsinf    = $1) if ($line =~m/parsimony-informative characters\s*=\s*(\d+)/i);
	($parsinf    = "All") if ($line =~m/All characters are parsimony-informative/);
        ($numtrees   = $1) if ($line =~m/trees retained\s*=\s*(\d+)/i);
        ($length     = $1) if ($line =~m/found\s*=\s*(\d+)/i);
	
	if ($line =~m/not enforced/i){
	    $conv = "none";
	}
	
	if ($line =~m/constraint-tree/){
	    if ($line =~m/NOT/){
		$conv = "without";
	    }
	    else {
		$conv = "with";
	    }
	}
	
        # look for end record marker                                              
	if ($line =~m/==>end_record/){
	    $bsdata{$constraint}{$conv} =
		[$parsinf, $numtrees, $length];
	}
    }
    close (BS);
    return (\%bsdata);
}

sub parse_ndi{
    my $self   = shift;
    my $file   = shift;

    my $partition;
    my $length;

    my $ndidata;
    open (NDI, "$file");

    while (my $line = <NDI>){
        chomp $line;

        # mine the data we need                                                   
        if ($line =~m/==>partition:\s*(\w+)/i){
            $partition = $1;
        }
	elsif ($line =~m/found\s*=\s*(\d+)/i){
	    $length = $1;
	}
        elsif ($line =~m/==>end_record/){
            $ndidata->{$partition} = $length;
	}
	elsif ($line =~m/Execution terminated due to errors/){
	    $ndidata = {};
	    last;
	}
        else {
            next;
        }
    }
    close (NDI);
    return ($ndidata);
}

sub parse_pbs{
    my $self   = shift;
    my $file   = shift;

    my $partition;
    my $length;
    my $total;
    my $elements;
    my @lengths;

    my $pbsdata;
    open (PBS, "$file");

    my $lencount = 0;
    while (my $line = <PBS>){
        chomp $line;

        # mine the data we need                                                                     
        if ($line =~m/==>partition:\s*(\w+)/i){
            $partition = $1;
        }
        elsif ($line =~m/^Length\s+/){
            my @lcl = split (/\s+/, $line);
            shift @lcl;
            push (@lengths, @lcl);
#            @lengths = split (/\s+/, $line);                                                                    
#           shift @lengths;                                                                                      
        }
        elsif ($line =~m/==>end_record/){
            next if ($partition eq "all");

            $total   += $_ foreach @lengths; # a bit iffy -- what if lengths go into next line?? possible??      
            $elements = @lengths;
            $length   = $total / $elements;

            $pbsdata->{$partition} = $length;
            @lengths  = ();
            $total    = 0;
            $elements = 0;
        }
        elsif ($line =~m/Execution terminated due to errors/){
            $pbsdata = {};
            last;
        }
        else {
            next;
        }
    }
    close (PBS);
    return ($pbsdata);
}

sub get_nchar{
    my $self = shift;
    return $self->{'alignment'}->{'nchar'};
}

sub get_charsets{
    my $self = shift;
    return $self->{'alignment'}->{'charset'};
}

sub get_root{
    my $self = shift;
    return $self->{'configuration'}->{'ROOT'};
}

sub get_config{
    my $self = shift;
    return $self->{'configuration'};
}

sub get_leaves{
    my $self = shift;
    
    my %nodes;
    my $lineages = $self->{'tree'};
    foreach my $node (keys %$lineages){
	foreach my $lineage (keys %{$lineages->{$node}}){
	    push @{$nodes{$node}}, @{ $lineages->{$node}->{$lineage} };
	}
    }

    return (\%nodes);
}

sub get_lineages{
    my $self = shift;
    return $self->{'tree'};
}

sub run_paup{
    my $self   = shift;
    my $outdir_cmds = shift;
    my $type   = shift;
    my $var    = shift;
    my $enf    = shift;
    
    if ($type eq "bs"){
	`$self->{'configuration'}->{PAUP} -n $outdir_cmds/$var.bs_commands.nex`;
    }
    elsif ($type eq "pbs"){
	`$self->{'configuration'}->{PAUP} -n $outdir_cmds/$var.$enf.pbs_commands.nex`;
    }
    elsif ($type eq "ndi"){
	`$self->{'configuration'}->{PAUP} -n $outdir_cmds/$var.$enf.ndi_commands.nex`;
    }
    else {
	print STDERR "No PAUP run: analysis type unknown\n";
    }
}

sub generate_partitions{
    my $self   = shift;
    my $window = shift;
    my $motion = shift;
    my $st     = shift;
    my $nchar  = shift;

    my %parts;
    my $counter = 0;
    my $i = $st;
    while (1){
	$counter++;
	my $start = $i;
	my $end   = $i + ($window - 1);

	my $coords = $start . "-" . $end;

	if ($end >= $nchar){
            $end = $nchar;
            $coords = $start . "-" . $end;
            $parts{$counter} = $coords;
            last;
	}

        $parts{$counter} = $coords;
	$i += $motion;
    }

    # over-ride previous charsets if defined in the matrix
    $self->{'alignment'}->{'charset'} = \%parts;

    return (\%parts);
}


1;
