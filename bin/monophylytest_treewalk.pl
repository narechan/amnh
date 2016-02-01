#!/usr/bin/perl -w

#####SETUP#####

# NOTE: needs at least bioperl 1.5.2

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;

my ($help, $treefile, $inlistfile, $outlistfile, $rootfile, $exceptions, $bootcut, $bootcollapse);
GetOptions(
	   'h|help'          => \$help,
	   't|treefile=s'   => \$treefile,
	   'a|outlistfile=s'       => \$outlistfile,
	   'b|inlistfile=s' => \$inlistfile,
	   'r|rootfile=s' => \$rootfile,
	   'e|exceptions=i' => \$exceptions,
	   'x|bootcut=i'    => \$bootcut,
	   'y|bootcollapse=i' => \$bootcollapse
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $inlistfile, $outlistfile, $exceptions, $bootcut, $bootcollapse){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

print STDERR "Working on $treefile\n";

# store the list of taxa outgroups
my $outlist = {};
open (L, "$outlistfile");
while (my $line = <L>){
    chomp $line;
    $outlist->{$line} = 1;
}
close (L);

# also store the ingroups                                                        
my $inlist = {};
open (I, "$inlistfile");
while (my $line = <I>){
    chomp $line;
    $inlist->{$line} = 1;
}
close (I);

# get the tree object and root if necessary
my $in_tree  = Bio::TreeIO->new(-file   => $treefile,
				-format => "newick",
				-internal_node_id => 'bootstrap');
my $tree     = $in_tree->next_tree;

# read the rootfile
my $roots = {};
open (R, "$rootfile");
while (my $line = <R>){
    chomp $line;
    my ($oid, $taxa, $contig) = split (/\t/, $line);
    my $root = $taxa . "#" . $contig;
    $roots->{$oid} = $root;
} 

my $treefilename;
if ($treefile =~/\//g){
    $treefile =~m/.*\/(.*)$/;
    $treefilename = $1;
}

else {
    $treefilename = $treefile;
}

my @treefile = split (/\./, $treefilename);
my $root = $roots->{$treefile[1]};
my @nodes = $tree->find_node(-id => $root);                                                                  
my $node = $nodes[0];                                                                                        
$tree->reroot($node);


#### OUTGROUP MINING ####

# find all outgroup clades
my @outseqskept;
my @outleaveskept;
my @outleaves = $tree->get_leaf_nodes;
foreach my $outleaf (@outleaves){
    my $id = $outleaf->id;
    my ($taxon, $acc) = split (/\#/, $id);
    if (exists ($outlist->{$taxon})){
        push (@outseqskept, $id);
        push (@outleaveskept, $outleaf);
    }
    else{
        next;
    }
}

# for every leaf walk up the tree until we hit a node                                                        
# where the opposite lineage contains something other than                                                     
# our taxa of interest and collapse the node                                                                                           
my $outclades = {};
foreach my $outleaf (@outleaveskept){
    my $leafid = $outleaf->id;

    # get the leaf's path all the down to the root                                                         
    my @path;
    while (defined $outleaf){
        push (@path, $outleaf);
        $outleaf = $outleaf->ancestor;
    }

    # take the leaf off the path                                                                            
    my $pathcache = {};
    my $bootcache = {};
    my $leafnode = shift @path;
    my $nodecounter = 1;
    $pathcache->{$nodecounter} = $leafnode;

    # query the descendants of each node up the path                                                    
    my $largestnode;
    foreach my $node (@path){
	$nodecounter++;
        my @lineages = $node->each_Descendent;

        # check to see if we're at the base, where there                                                    
        # is only one lineage (discounts the outgroup)                                                       
        # outgroup will never be included in a clade                                                         
        my $lincount = @lineages;
        if ($lincount <= 1){
	    my $previousnode = $nodecounter - 1;
	    $largestnode = $pathcache->{$previousnode};
            last;
        }

	# skip the lineage that has the leaf of interest                                                       
        my $keptlineage;
        foreach my $lineage (@lineages){
            my @children;

            if ($lineage->is_Leaf){
                push (@children, $lineage);
            }
            else {
                my @descendents = $lineage->get_all_Descendents;
                @children = grep{$_->is_Leaf} @descendents;
            }

            my $insignal = 0;
            foreach my $child (@children){
                my $id = $child->id;
                if ($id eq $leafid){
                    $insignal++;
                    last;
                }
                else{
                    next;
                }
            }
            if ($insignal > 0){
                next;
            }
            else{
                $keptlineage = $lineage;
                last;
            }
        }

	# get the children in the kept lineage                                                                
        my @children;
        if ($keptlineage->is_Leaf){
            push (@children, $keptlineage);
        }
        else {
            my @descendents = $keptlineage->get_all_Descendents;
            @children = grep{$_->is_Leaf} @descendents;
        }

        # log all ingroups in the kept lineage                                                               
        my @ins;
        foreach my $leaf (@children){
            my $id = $leaf->id;

            my ($taxon, $acc) = split (/\#/, $id);
            if (exists($inlist->{$taxon})){
                push (@ins, $taxon);
            }
            else {
                next;
            }
        }

	# find the bootstrap
	my $bootnode = $node->bootstrap;
	
	# if ingroup encountered, log the pathcache as complete
	# and search for the largest outgroup clade along the 
	# path that satisfies the bootstrap collapse criterion
	my $ins = @ins;
	if ($ins > 0){

	    # check to see if the leaf is the only thing in there;                   
            # if so, assign it as the largest node and skip the boot test            
            my @pathcache = keys %$pathcache;
            my $pathcachecnt = @pathcache;
            if ($pathcachecnt == 1){
                $largestnode = $pathcache->{1};
            }
	    else{
		foreach my $nodenum (sort {$b <=> $a} keys %$pathcache){
                    # we've drilled all the way down to the leaf,                 
		    # so take that as the largest node                          
		    
                    if ($nodenum == 1){
                        $largestnode = $pathcache->{$nodenum};
                        last;
                    }
                    elsif ($bootcache->{$nodenum} >= $bootcollapse){
                        $largestnode = $pathcache->{$nodenum};
                        last;
                    }
                    else{
                        next;
                    }
                }
            }
            last;
        }
	else {
	    my $bootstrap = $node->bootstrap;
            $pathcache->{$nodecounter} = $node;
            $bootcache->{$nodecounter} = $bootstrap;
        }
	
    }

    # get all the taxa in the largest node and record                                                         
    my @largedescendents = $largestnode->get_all_Descendents;
    my @largechildren = grep{$_->is_Leaf} @largedescendents;
    my @llseqs;

    # get the bootstrap (in this case, coded as node id)                                                      
    my $bootstrap = $largestnode->bootstrap;
    
    # if there are no descendants, this is a singleton clade                                              
    # otherwise print all leaves to define the clade                                                         
    # unique-ify the species in each collapsed outgroup clade
    my $numchildren = @largechildren;
    if ($numchildren == 0){
	my ($sp, $acc) = split (/\#/, $leafid);
	$outclades->{$leafid} = $leafid;
    }
    else{
	foreach my $largeleaf (@largechildren){                                                               
	    my $llid = $largeleaf->id;                                                                        
	    push (@llseqs, $llid);                                                                            
	    my ($sp, $acc) = split (/\#/, $llid);                                                             
	}                                                                                                        
	my $llseqs = join ",", sort @llseqs;                                                                  
	$outclades->{$leafid} = $llseqs; 
    }
}

#### INGROUP MINING ####

# get all the leaf nodes and screen for the ones
# with our taxa of interest
my @seqskept;
my @leaveskept;
my @leaves = $tree->get_leaf_nodes;
foreach my $leaf (@leaves){
    my $id = $leaf->id;
    my ($taxon, $acc) = split (/\#/, $id);
    if (exists ($inlist->{$taxon})){
	push (@seqskept, $id);
	push (@leaveskept, $leaf);
    }
    else{
	next;
    }
}

# for every leaf walk up the tree until we hit a node
# where the opposite lineage contains something other than 
# our taxa of interest
my $clades = {};
my $monoclades = {};
foreach my $leaf (@leaveskept){
    my $leafid = $leaf->id;

    # get the leaf's path all the down to the root
    my @path;
    while (defined $leaf){
	push (@path, $leaf);
	$leaf = $leaf->ancestor;
    }

    # take the leaf off the path
    my $pathcache = {};
    my $bootcache = {};
    my $leafnode = shift @path;
    my $nodecounter = 1;
    $pathcache->{$nodecounter} = $leafnode;

    # query the descendants of each node up the path
    my $largestnode;
    my $largestmononode;
    my $mononodelock = 0;
    my $outseen = {};
#    my @inspath;
    my @monoinspath;
    my $leafboot;
    foreach my $node (@path){
	$nodecounter++;
	($leafboot = $node->bootstrap) if ($nodecounter == 2);
	my @lineages = $node->each_Descendent;

	# check to see if we're at the base, where there 
	# is only one lineage (discounts the outgroup)
	# outgroup will never be included in a clade
	my $lincount = @lineages;
	if ($lincount <= 1){
	    my $previousnode = $nodecounter - 1;
            $largestnode = $pathcache->{$previousnode};
	    last;
	}

	# skip the lineage that has the leaf of interest
	my $keptlineage;
	foreach my $lineage (@lineages){
	    my @children;

            if ($lineage->is_Leaf){
                push (@children, $lineage);
            }
            else {
                my @descendents = $lineage->get_all_Descendents;
                @children = grep{$_->is_Leaf} @descendents;
            }

	    my $insignal = 0;
	    foreach my $child (@children){
		my $id = $child->id;
		if ($id eq $leafid){
		    $insignal++;
		    last;
		}
		else{
		    next;
		}
	    }
	    if ($insignal > 0){
		next;
	    }
	    else{
		$keptlineage = $lineage;
		last;
	    }
	}
	    
	# for the opposite lineage, check if there are any outgroup taxa,
	# if so, find the outgroup's clade based on outgroup mining and log.
	# if more than x clade exceptions found, bail.
	
	# get the children in the kept lineage
	my @children;
	if ($keptlineage->is_Leaf){
	    push (@children, $keptlineage);
	}
	else {
	    my @descendents = $keptlineage->get_all_Descendents;
	    @children = grep{$_->is_Leaf} @descendents;
	}

	# log all outgroups in the kept lineage and get their unique-ified clades
	# log all ingroups for the ingroups stream
	my @outs;
	my $ins;
	my $outstokeep;
	foreach my $leaf (@children){
            my $id = $leaf->id;
            my ($taxon, $acc) = split (/\#/, $id);
	    if (exists($outlist->{$taxon})){
		push (@outs, $taxon);
		$outstokeep->{$outclades->{$id}} = 1;
	    }		
	    elsif (exists($inlist->{$taxon})){
		$ins .= $taxon;
	    }
	    else {
		next;
	    }
        }
#	push (@inspath, $ins);
	(push (@monoinspath, $ins)) unless ($mononodelock > 0);

	# do monophyly specific stuff (stuff to do when first outgroup is seen)
	my $outs = @outs;
	if (($outs > 0) and ($mononodelock == 0)){
	    $mononodelock++;

	    pop @monoinspath;
	    my $monoinspath = join ":", @monoinspath;
	    if ($leafboot){
		print STDERR "$treefilename\t$leafid\t$leafboot\t$monoinspath\n";
	    }
	    else {
		print STDERR "$treefilename\t$leafid\tNA\t$monoinspath\n";
	    }
	    my @pathcache = keys %$pathcache;
	    my $pathcachecnt = @pathcache;
	    if ($pathcachecnt == 1){
		$largestmononode = $pathcache->{1};
	    }
	    else {
		foreach my $nodenum (sort {$b <=> $a} keys %$pathcache){
		    # we've drilled all the way down to the leaf,                                             
		    # so take that as the largest node                                                        
		    if ($nodenum == 1){
			$largestmononode = $pathcache->{$nodenum};
			last;
		    }
		    elsif ($bootcache->{$nodenum} >= $bootcut){
			$largestmononode = $pathcache->{$nodenum};
			last;
		    }
		    else{
			next;
		    }
		}
	    }
	}

	# collate all the outgroups seen to this point, collpasing species 
	# in outgroup clades
	foreach my $clade (keys %$outstokeep){
	    my @clade = split (/\,/, $clade);
	    my $spuniq = {};
	    foreach my $species (@clade){
		my ($taxon, $acc) = split (/\#/, $species);
		$spuniq->{$taxon}++;
	    }
	    foreach my $species (keys %$spuniq){
		$outseen->{$species}++;
	    }
	}
	
	# if any outgroup species seen more than the allowable 
	# number of exceptions, then bail
	my $sig = 0;
	foreach my $spout (keys %$outseen){
	    if ($outseen->{$spout} > $exceptions){
		$sig++;
		last;
	    }
	}
	
	# if we've hit a exceptions violation, do not log the node, but
	# go back down the node tree until we find a node that statsifies 
	# the bootstrap criterion
	if ($sig == 1){

	    # check to see if the leaf is the only thing in there;
	    # if so, assign it as the largest node and skip the boot test
	    my @pathcache = keys %$pathcache;
	    my $pathcachecnt = @pathcache;
	    if ($pathcachecnt == 1){
		$largestnode = $pathcache->{1};
	    }
	    else {
		foreach my $nodenum (sort {$b <=> $a} keys %$pathcache){
		    # we've drilled all the way down to the leaf,                                           
                    # so take that as the largest node                                                        
                    if ($nodenum == 1){
			$largestnode = $pathcache->{$nodenum};
                        last;
		    }
		    elsif ($bootcache->{$nodenum} >= $bootcut){
			$largestnode = $pathcache->{$nodenum};
			last;
		    }
		    else{
			next;
		    }
		}
	    }
	    last;
	}
	else {
	    my $bootstrap = $node->bootstrap;
	    $pathcache->{$nodecounter} = $node;
	    $bootcache->{$nodecounter} = $bootstrap;
        }
    }
    
    # get the exceptions path
#    pop @inspath;
#    my $inspath = join ":", @inspath;
#    print STDERR "$leafid\t$inspath\n";
    
    ### exceptions monophyly harvesting ###
    
    # get all the taxa in the largest node and record                                                       
    my @largedescendents = $largestnode->get_all_Descendents;
    my @largechildren = grep{$_->is_Leaf} @largedescendents;
    my @llseqs;
    
    # get the bootstrap (in this case, coded as node id)
    my $bootstrap = $largestnode->bootstrap;

    # if there are no descendants, this is a singleton clade                                              
    # otherwise print all leaves to define the clade                                                        
    my $numchildren = @largechildren;
    if ($numchildren == 0){
        $clades->{$leafid} = "NA";
    }
    else{
        foreach my $largeleaf (@largechildren){
            my $llid = $largeleaf->id;
            push (@llseqs, $llid);
        }
        my $llseqs = join ",", sort @llseqs;
        $clades->{$llseqs} = $bootstrap;
    }
    
    ### strict monophyly harvesting ###

    # deal with the case where the root is the only outgroup
    unless (%$outseen){
	$largestmononode = $largestnode;
    }
    
    # get all the taxa in the largest strict mono node and record                                         
    my @largemonodescendents = $largestmononode->get_all_Descendents;
    my @largemonochildren = grep{$_->is_Leaf} @largemonodescendents;
    my @monollseqs;

    # get the bootstrap (in this case, coded as node id)                                            
    my $monobootstrap = $largestmononode->bootstrap;

    # if there are no descendants, this is a singleton clade                                                 
    # otherwise print all leaves to define the clade                                                       
    my $nummonochildren = @largemonochildren;
    if ($nummonochildren == 0){
        $monoclades->{$leafid} = "NA";
    }
    else{
        foreach my $largemonoleaf (@largemonochildren){
            my $monollid = $largemonoleaf->id;
            push (@monollseqs, $monollid);
        }
        my $monollseqs = join ",", sort @monollseqs;
        $monoclades->{$monollseqs} = $monobootstrap;
    }
}

### exceptions monophyly printing ###

# print all the unique largest clades that                                                              
# obey the largest monophyly rules                                                                       
my $cladecounter = 0;
foreach my $clade (keys %$clades){
    $cladecounter++;
    if ($clades->{$clade}){
	print "$treefilename\tEXCEP\tclade$cladecounter\t$clade\t$clades->{$clade}\n";
    }
    else {
	print "$treefilename\tEXCEP\tclade$cladecounter\t$clade\tNULL\n";
    }
}

### strict monophyly printing (to STDERR) ###

# print all the unique largest clades that                                                     
# obey the largest monophyly rules                                                                  
my $monocladecounter = 0;
foreach my $monoclade (keys %$monoclades){
    $monocladecounter++;
    if ($monoclades->{$monoclade}){
        print "$treefilename\tMONO\tclade$monocladecounter\t$monoclade\t$monoclades->{$monoclade}\n";
    }
    else {
        print "$treefilename\tMONO\tclade$monocladecounter\t$monoclade\tNULL\n";
    }
}
























### OLD CODE that does more sophisticated things like finding the LCA in the 
# query clade to check for monophyly of outgroups. monophyly may only be working 
# for the case where there is only 1 exception.

=head

            # for the opposite lineage check if the taxa are all ingroup,                                             
        # all outgroup, or a mixture. if all ingroup, accept and move up                                          
        # a node. if all outgroup, record exeption, accept and move up a node.                                    
        # if a mixture, drill and check whether the outgroup taxa are                                             
        # monophyletic. if so, record exception, accept and monve up a node.                                      
        # if we've encountered our second outgroup regardless of monophyly,                                       
        # bail and record                                                                                         


	# ennumerate the species
	my $species;
	my @out;
	foreach my $leaf (@children){
	    my $id = $leaf->id;
	    my ($taxon, $acc) = split (/\#/, $id);
	    $species->{$taxon}++;
	    (push (@out, $leaf)) if (exists($outlist->{$taxon}));
	}
	
	# check ingroup and outgroup status
	my $in = 0;
	my $out = 0;
	foreach my $sp (keys %$species){
	    if (exists($inlist->{$sp})){
		$in++;
	    }
	    else {
		$out++;
	    }
	}

	# check only for outgroup scenarios. if all in lineage
	# is in the ingroup, update patchcache and go to the 
	# next node (so do nothing)

	# if lineage is all ingroup except for 1 outgroup,
	# increase signal by 1
	if (($in >= 1) and ($out == 1)){
	    $signal++;
	}

	# if lineage is all outgroup (by def monophy), increase signal by 1
	if (($in == 0) and ($out >= 1)){
	    $signal++;
	}

	# if lineage is mixed, and more than one outgroup check to see if outgroup 
	# lineage is monophyletic
	if (($in >= 1) and ($out > 1)){
	    $signal++;
	    my $lca = $tree->get_lca(-nodes => \@out);

	    my @lcachildren;
	    if ($lca->is_Leaf){
		push (@lcachildren, $lca);
	    }
	    else {
		my @lcdescendents = $lca->get_all_Descendents;
		@lcachildren = grep{$_->is_Leaf} @lcdescendents;
	    }

	    my $lcspecies = {};
	    foreach my $lc (@lcachildren){
		my $id = $lc->id;
		my ($taxon, $acc) = split (/\#/, $id);
		$lcspecies->{$taxon}++;
	    }
	    
	    my $lsignal = 0;
	    foreach my $lsp (keys %$lcspecies){
		if (exists($inlist->{$lsp})){
		    $lsignal++;
		}
		else {
		    next;
		}
	    }

	    # if not monophyletic, cause an immediate bail and report
	    # NOTE: this may only be appropriate for the one exception case!!
	    if ($lsignal > 0){
		$signal = 100;
	    }
	}
	
	if ($signal >= ($exceptions + 1)){
            $largestnode = $pathcache;
            last;
        }
        else {
            $pathcache = $node;
        }
    }

    # get all the taxa in the largest node and record                                                      
    my @largedescendents = $largestnode->get_all_Descendents;
    my @largechildren = grep{$_->is_Leaf} @largedescendents;
    my @llseqs;

    # if there are no descendants, this is a singleton clade                                                 
    # otherwise print all leaves to define the clade                                                        
    my $numchildren = @largechildren;
    if ($numchildren == 0){
        $clades->{$leafid} = 1;
#       print "$leafid\t$leafid\n";                                                                        
    }
    else{
        foreach my $largeleaf (@largechildren){
            my $llid = $largeleaf->id;
            push (@llseqs, $llid);
        }
        my $llseqs = join ",", sort @llseqs;
        $clades->{$llseqs} = 1;
#       print "$leafid\t$llseqs\n";                                                                        
    }
}

# print all the unique largest clades that                                                              
# obey the largest monophyly rules                                                                          
my $cladecounter = 0;
foreach my $clade (keys %$clades){
    $cladecounter++;
    print "clade$cladecounter\t$clade\n";
}

=cut
