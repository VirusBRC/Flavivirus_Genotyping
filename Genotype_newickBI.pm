package Genotype_newickBI;

# This file takes a newick phylogeny tree, and the name of a node, find if the node is located with a branch of the tree of same genotype
# Requires access to:
#
#use strict;
use warnings;
use File::Spec::Functions;
use IO::String;
use Bio::SeqIO;
use POSIX qw(strftime ceil floor);
use Bio::TreeIO;
use Bio::AlignIO;
#use Bio::SimpleAlign;
use Scalar::Util 'refaddr';
use Bio::Tools::Run::StandAloneBlast;

use English;
use Carp;
use Data::Dumper;
use version; our $VERSION = qv('1.3.0'); # Jun 06, 2013
#use Getopt::Long;
use Annotate_Def;
use Genotype_util;

my $debug_all = Genotype_util::getDebugAll();

############## Constants ################
# To add a new species:
# 0. get all refseqs for the genotype in fasta file, make sure defline is >1a.AB123456
# 1. create refseq MSA: muscle -in refseq_stlouisenvephalitis_33.fasta -out refseq_stlouisenvephalitis_33.aln -stable
# 2. add the $RMSA_species to following 2 lists: $RMSA_<species> and $RMSA hash
# 3. add any refseq to Genotype_refseq_def.txt, and define the shorthand
# 4. add any long refseq to refseq_Flaviviridae<nn>.faa, and change the file name.
#    To create database needed by blast, run: formatdb -i refseq_Flaviviridae84.faa -p F
# Note: The automatic detection of species is for convenience, but not a fool proof method, user input is preferred.
# The St Louis and Bovine diarrhea are especially difficult.


#    my $newHCV_ICTV = 1; # Sep 30, 2013
#    my $fmtGenotype = '\d[-._\dA-Za-z]*';
#    $fmtGenotype = '([0-9]\w+|\d-[A-Za-z]+|NA|.+)' if (!$newHCV_ICTV);
    my $fmtGenotype = Genotype_Def::getfmtGenotype;

############## Subroutines ################

=head2 newickBI
Calculate the Branching Index (BI) of a target node from a tree file
# This code was first developed at LANL, then adapted by VBRC
# At VIPR, improvements include handling of intermediate nodes, as well as re-rooting the tree
# to handle cases where the target is very close to tree root.
=cut

sub newickBI {
    my ($treeFile, $focalNode, $taxon, $debug) = @_;

    $debug = 0 || $debug_all;
    my $subname = 'Genotype_newickBI::newickBI';
    my $dumper_debug = 0;

$debug && print STDERR "$subname: \$treeFile=$treeFile\n";
$debug && print STDERR "$subname: Tree file=$treeFile\n";
$debug && print STDERR "$subname: focalNode is '$focalNode'\n";
#    if ($focalNode =~ m/_$/) {
#        $focalNode =~ s/_$/./i; 
#        $debug && print STDERR "$subname: \$focalNode changed to: '$focalNode'\n";
#    }
die "$subname: ERROR: cannot find $treeFile ($!)" if !-e $treeFile;
die "$subname: ERROR: please specify a focal taxon" if (!$focalNode);

my $input = new Bio::TreeIO(-file => $treeFile, -format => "newick");
#my $input = new Bio::TreeIO(-file => $treeFile, -format => "nexus"); # "newick";

my $PRESERVE_FAILURES = 1; # log errors for forensic analysis
my $OMIT_ORPHANS = 1;  # requires taxon branch to be at subtype level not above
my $PRINT_TAXON_ID = 1;
my $myTaxonID = "NA";
my $failDir = "./nbi-failures" if $PRESERVE_FAILURES;

my $newickBI_result = '';
my $itree = -1;
while (my $tree = $input->next_tree) {
    $itree++;
    $dumper_debug && print STDERR "$subname: \$itree=$itree \$tree='\n". Dumper($tree) . "End of \$tree\n\n";
    my @leafs = $tree->get_leaf_nodes;
    my @nodes = $tree->get_nodes;
    $debug && print STDERR "$subname: \$itree=$itree has \$#nodes=$#nodes\n";
    next if ($#nodes<=0);

    my $BI = 0;
    my $myTaxonID = "NA";
    my $sub = "NA";

    my $focus = Genotype_newickBI::find_focalNode( $tree, $focalNode);
    # Re-root the tree
        Genotype_newickBI::re_root_tree( $tree, $focus, $debug);
        $focus = Genotype_newickBI::find_focalNode( $tree, $focalNode); # find the focalNode again after re_root
    $debug && print STDERR "$subname: Focus is: " . $focus->to_string . "\n";

    my $parent = $focus->ancestor;
    $dumper_debug && print STDERR "newickBI:\$parent=".Dumper($parent)."\n";
#    &mydie("ERROR: $focalNode has no parent w/ $#leafs leafs") if !$parent;
    # BI is defined as zero if there is no parent
    if (!$parent) {
        $debug && print STDERR "$subname: \$BI=0 because 'there is no parent'\n";
        $BI=0;
        $myTaxonID ="NA";
        my $msg = sprintf "BI=%5.3f", $BI;
        $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
        $msg .= sprintf "\t$focalNode has no parent w/ $#leafs leafs\n";
        $newickBI_result .= $msg;
        next;
    }


    my $grandparent = $parent->ancestor;
    $dumper_debug && print STDERR "newickBI:\$grandparent=".Dumper($grandparent)."\n";
#    &mydie("ERROR: $focalNode has no grandparent w/ $#leafs leafs") if !$grandparent;
    # BI is defined as zero if there is no grandparent
    if (!$grandparent) {
        $debug && print STDERR "$subname: \$BI=0 because 'there is no grandparent'\n";
        $BI=0;
        $myTaxonID ="NA";
        my $msg = sprintf "BI=%5.3f", $BI;
        $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
        $msg .= sprintf "\t$focalNode has no grandparent w/ $#leafs leafs\n";
        $newickBI_result .= $msg;
        next;
    }


    # Get all descendents, including all end-nodes, and nodes in the middle
    my @sibNodes = $parent->get_all_Descendents;
    $debug && print STDERR "$subname: \$parent \@sibNodes=$#sibNodes\n";
#    for my $isib (0 .. $#sibNodes) {
#        my $sib = $sibNodes[$isib];
#        $dumper_debug && print STDERR "$subname: \$isib=$isib \$sib='\n". Dumper($sib) . "End of \$sib\n\n";
#    }
    my $nSibs=0;
    my $sister;
#    my $BI=0;
    my @sibIDs = (); # The genotype from each sibling
    my @extIDs = (); # The genotype from node that is not a sibling
    my $distanceToNearestSibling = 99999999;

    foreach my $isib (0 .. $#sibNodes) {
        my $sib = $sibNodes[$isib];
        $debug && print STDERR "$subname: \$isib=$isib \$focus=".refaddr($focus)." \$sib=".refaddr($sib)."\n";
        next if ($sib eq $focus);
#        $dumper_debug && print STDERR "\$isib=$isib \$sib=\n" . Dumper($sib) . "\n";

        if (defined ($sib->id) && $sib->id) {
            $debug && print STDERR "$subname: \$isib=$isib id=".$sib->id." \$focus=".refaddr($focus)." \$sib=".refaddr($sib)."\n";
            if ($sib->id =~ m/^\s*0[.]000000\s*$/) {
                $debug && print STDERR "$subname: \$isib=$isib id=".$sib->id." Skipped.\n";
                next;
            }
            $debug && print STDERR "$subname: \$taxon=$taxon \$isib=$isib \$sib->id='". $sib->id ."'\n";

#            if ($sib->id =~ /^[0-9.]+$/ || $sib->id =~ /^\s*$/) {
#                $debug && print STDERR "$subname: \$isib=$isib Found sibling with id as all numbers. Skip.\n";
##                next; # trying to handle a sibling with id=0.976000
#            }
#            my @currentClade = split /\./, $sib->id;
            my @currentClade = split /\|/, $sib->id;
#            my $genotype_pattn = '^(\d[a-z]*|\d-[A-Za-z]+)$';
            my $genotype_pattn = '^\d[-\d._A-Za-z]*$';
            if ($taxon =~ /^NoV/i) {
                @currentClade = split /[|]/, $sib->id;
                $genotype_pattn = '^\w*[.]\w+$';
            }
            $debug && print STDERR "$subname: \$isib=$isib \@currentClade='@currentClade' \$genotype_pattn='$genotype_pattn'\n";
            $currentClade[0] =~ s/'//;         #' restore syntactical highlighting
#            push @sibIDs, $currentClade[0];
            if ($currentClade[0] =~ m/($genotype_pattn)/i) {
                my $seen = 0;
                for my $sib (@sibIDs) { $seen = 1 if ($sib eq $currentClade[0]); }
                push @sibIDs, $currentClade[0] if (!$seen);
                $debug && print STDERR "$subname: \$1=$1 \$isib=$isib \@sibIDs='@sibIDs'\n";
            } else {
                $debug && print STDERR "$subname: Problem: \$isib=$isib \$currentClade[0]='$currentClade[0]'\n";
                next; # skip if there is a problem
            }

            # find nearest sibling and store its subtype
            my $distanceToSibling = $tree->distance(-nodes => [$focus,$sib]);
            $debug && print STDERR "$subname: \$isib=$isib \$distanceToSibling=$distanceToSibling\n";

            if ($distanceToSibling < $distanceToNearestSibling) {
                if ($currentClade[0] =~ m/$genotype_pattn/i) {
                    $distanceToNearestSibling = $distanceToSibling;
                    $myTaxonID = $currentClade[0];
                    $sub = $currentClade[1] if ($taxon =~ /^NoV/i && $myTaxonID =~ m/^II.p*4$/i);
                }
            }
            undef @currentClade;
            $debug && print STDERR "$subname: \$distanceToSibling=$distanceToSibling\n";
            $debug && print STDERR "$subname: \$distanceToNearestSibling=$distanceToNearestSibling\n";
            $debug && print STDERR "$subname: \$myTaxonID=$myTaxonID\n";
        }

        $dumper_debug && print STDERR "$subname:\$sib->ancestor=".Dumper($sib->ancestor)."\n";
        $dumper_debug && print STDERR "$subname:\$parent=".Dumper($parent)."\n";
        $debug && print STDERR "$subname: \$sib->ancestor=".refaddr($sib->ancestor)." id='".(($sib->ancestor->id) ? $sib->ancestor->id : '')."'\n";
        $debug && print STDERR "$subname:        \$parent=".refaddr($parent)." id='".(($parent->id) ? $parent->id : '')."'\n";
        if ($sib->ancestor eq $parent) {
            $sister=$sib;
            $nSibs++;
        }
        $debug && print STDERR "$subname: \$nSibs=$nSibs \@sibIDs='@sibIDs'\n";
    }

    $debug && print STDERR "$subname: \@sibIDs='@sibIDs'\n\n";
    my $areSibsOneSubtype = 1;
if ( $taxon ne 'ZIKA' ) { # Original code
    foreach my $sibA (@sibIDs) {
    # @sibIDs='1b 1b'
        foreach my $sibB (@sibIDs) {
            $debug && print STDERR "$subname: \$sibA='$sibA' \$sibB='$sibB'\n";
            $areSibsOneSubtype = 0 if ($sibA ne $sibB);
            last unless $areSibsOneSubtype;
        }
        last unless $areSibsOneSubtype;
    }
} else { # try to roll up the genotype so that 1a,1b => 1; but this doesn't work -02/05/2016
    foreach my $isibA (0 .. $#sibIDs) {
        my $sibA = $sibIDs[$isibA];
    # @sibIDs='1b 1b'
        foreach my $isibB (0 .. $#sibIDs) {
            my $sibB = $sibIDs[$isibB];
            $debug && print STDERR "$subname: \$sibA='$sibA' \$sibB='$sibB'\n";
            next if ($sibA eq $sibB);
            if ($sibA =~ m/^(\d+)[a-z]+/i) {
              my $p = $1;
              if ($sibB =~ m/^($p)[a-z]*/) {
                $sibIDs[$isibA] = $p; # changes $sibA to number only, such as 1 in 1a -02/05/2016
                $sibIDs[$isibB] = $p; # changes $sibA to number only, such as 1 in 1a -02/05/2016
              } else {
                $areSibsOneSubtype = 0;
              }
            }
            last unless $areSibsOneSubtype;
        }
        last unless $areSibsOneSubtype;
    }
    # Change the assignment, since all sibIDs should be same now if $areSibsOneSubtype=1
    $myTaxonID = $sibIDs[0] if ($areSibsOneSubtype && $#sibIDs>=0);
}
    $debug && print STDERR "$subname: \$myTaxonID=$myTaxonID \@sibIDs='@sibIDs'\n\n";
    $debug && print STDERR "$subname: \$areSibsOneSubtype=$areSibsOneSubtype\n\n";

    # BI is defined as zero if we are clearly more basal than one subtype-specific clade
    # i.e., if not all the siblings belong to same genotype, such as 1a 
    # Or if there is no element in @sibIDs, such as in new Zika case
    if (!$areSibsOneSubtype || $#sibIDs<0) {
        $debug && print STDERR "$subname: \$BI=0 because 'clearly more basal than one subtype-specific clade, or \@sibIDs is empty'\n";
        $BI=0;
        $myTaxonID ="NA";
        my $msg = sprintf "BI=%5.3f", $BI;
        $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
        my $count = ($#sibIDs<=2) ? $#sibIDs : 2;
        my $delimiter = ',';
        $msg .= "\tsibIDs=". join($delimiter, @sibIDs[0 .. $count]);
        $msg .= '...' if ($#sibIDs>$count);
        $msg .= "\n";
        $newickBI_result .= $msg;
        $debug && print STDERR "$subname: \$taxon=$taxon result=$newickBI_result";
        next;
    }

    # locate "extended family" node/s that are neither sibs nor focus
    $dumper_debug && print STDERR "$subname:\$grandparent=".Dumper($grandparent)."\n";
    my @extNodes = $grandparent->get_all_Descendents;
    foreach my $iext (0 .. $#extNodes) {
        my $ext = $extNodes[$iext];
        $debug && print STDERR  "$subname: The extended family: \$iext=$iext '" . (($ext->id) ? $ext->id : '') ."'\n";
        next if $ext eq $focus;

        my $isExtNodeSibling = 0;
        foreach my $sib (@sibNodes) {
#            $debug && print STDERR  "\$sib=$sib \$focus=$focus\n";
            next if $sib eq $focus;
            next if (!($ext->id) || !($sib->id));
            $isExtNodeSibling = 1 if ($sib eq $ext);
        }
        $debug && print STDERR  "$subname: \$ext='".(($ext->id) ? $ext->id : '')."' \$isExtNodeSibling=$isExtNodeSibling\n";
        next if $isExtNodeSibling;

        # get the clade designation/label for extended family
        if (defined ($ext->id) && $ext->id) {
            $debug && printf STDERR "$subname: clade ext %s\n", $ext->id;

#            my @currentClade = split /\./, $ext->id;
            my @currentClade = split /\|/, $ext->id; # genotype identifier is separated by '|' - Oct 09, 2013
            $debug && print STDERR  "$subname: \@currentClade='@currentClade'\n";
            $currentClade[0] =~ s/'//;         #' restore syntactical highlighting
#            push @extIDs, $currentClade[0] if ($currentClade[0] =~ m/^(\d[a-z]*|\d-[A-Za-z]+)$/i);
            if ($currentClade[0] =~ m/^$fmtGenotype$/i) {
                push @extIDs, $currentClade[0];
                $debug && print STDERR  "$subname: \$currentClade[0]='$currentClade[0]' Found\n";
            } else {
                $debug && print STDERR  "$subname: \$currentClade[0]='$currentClade[0]' Problem\n";
            }
            undef @currentClade;
        }
    }
    $debug && print STDERR  "$subname: \@extIDs='@extIDs'\n\n";

    my $isSubtypeClade = 1; # set this to false in evidence to the contrary
    foreach my $extID (@extIDs) {
        # check if these are all the same, BI=1
        $debug && printf STDERR "$subname: \$extID='$extID'\n";

        foreach my $sibID (@sibIDs) {
            $debug && printf STDERR "$subname: extID=$extID sibID=$sibID\n";
            $isSubtypeClade = 0 if ($sibID ne $extID);
            $debug && printf STDERR "$subname: \$isSubtypeClade=$isSubtypeClade\n";
            last if !$isSubtypeClade;
        }
        last if !$isSubtypeClade;

#        $myTaxonID = $extID;
    }

# BI is defined as unity if we are clearly within a subtype-specific clade
    if ($isSubtypeClade) {
        $debug && print STDERR "$subname: \$BI=1 because 'clearly within a subtype-specific clade'\n";
        $BI=1;
        my $msg = sprintf "BI=%5.3f", $BI;
        $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
        $msg .= sprintf "\n";
        $newickBI_result .= $msg;
        next;
    }

    $debug && print STDERR "$subname: \$nSibs=$nSibs\n";
    $debug && print STDERR "$subname: \$sister=$sister\n";
    &mydie("ERROR: node $focalNode parent is polytomy so BI is undefined") if $nSibs>1;
    &mydie("ERROR: node $focalNode is an only child so BI is undefined") if $nSibs==0 or !$sister;

    my $alpha = $tree->distance(-nodes => [$grandparent,$parent]);
    ($alpha<0) && print STDERR "$subname: Warning: negative branch length in tree: \$alpha=$alpha\n";
    my $beta = $tree->distance(-nodes => [$sister,$parent]);
    ($beta<0) && print STDERR "$subname: Warning: negative branch length in tree: \$beta=$beta\n";

    # store taxon ID from sister
#    if (defined ($sister->id) && $sister->id) {
#	my @sisterClade = split /\./, $sister->id;
#	$sisterClade[0] =~ s/'//;         #' restore syntactical highlighting
#	$myTaxonID = $sisterClade[0];
#	undef @sisterClade;
#    }

    $debug && print STDERR "$subname: alpha=$alpha beta=$beta\n";

#    &mydie("ERROR: cannot divide by zero (a=$alpha, b=$beta myTaxonID=$myTaxonID)") if ($alpha+$beta == 0 && !$BI);
    if ($alpha+$beta == 0 && !$BI) {
        $BI = 0;
        $myTaxonID = 'NA';
        my $msg = sprintf "BI=%5.3f", $BI;
        $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
        $msg .= sprintf "\tZero distance between sibling and grandparent: a+b=0\n";
        $newickBI_result .= $msg;
        next;
    } else {
        $BI = $alpha/($alpha+$beta) unless $BI;
        $debug && print STDERR "$subname: \$alpha=$alpha \$beta=$beta \$BI=$BI\n";
    }
    $BI = 0 if $BI < 0;
    $BI = 1 if $BI > 1;

    $debug && print STDERR "$subname: \$PRINT_TAXON_ID=$PRINT_TAXON_ID\n";
    $debug && print STDERR "$subname: \$myTaxonID=$myTaxonID\n";
    my $msg = sprintf "BI=%5.3f", $BI;
    $msg .= sprintf "\tgenotype=%s", $myTaxonID if $PRINT_TAXON_ID;
        $msg .= "\tsub=$sub";
    $msg .= sprintf "\n";
    $newickBI_result .= $msg;
}

    $newickBI_result = "BI=\tgenotype=NA\tsub=$sub\n" if (!$newickBI_result);
    $debug && print STDERR "$subname: \$taxon=$taxon result='$newickBI_result'\n";
    $debug && print STDERR "$subname: leaving subroutine\n";
    return $newickBI_result;
} # sub newickBI


=head2 re_root_tree

Find an intermediate node that best separates major branches and re-root the tree to this node
 This is done by: find the lest populous groups, re-root with its parent node, then find
 the most populous groups, re-root with its parent node

=cut

sub re_root_tree {
    my ($tree, $focus, $debug) = @_;

    $debug = 0 || $debug_all;
    my $subname = 'Genotype_newickBI::re_root_tree';

    # Re-root the tree, using the parent node that covers different genotypes as root node
    my @leafs = $tree->get_leaf_nodes;
    $debug && print STDERR "$subname: \$focus=".$focus->id." \$#leafs=$#leafs\n";
    return if ($#leafs<2);
    my $new_root = '';
    if ( 1 ) {
        # First, determine the groups
        my $groups;
        $debug && print STDERR "$subname:";
        for my $il (0 .. $#leafs) {
            my $leaf_id = $leafs[$il]->id;
            $debug && print STDERR "\t$il=$leaf_id";
            push @{$groups->{$1}}, $leafs[$il] if ($leaf_id =~ /^[']*(\d|I+[.])/); # Sometimes ids have \'
#            push @{$groups->{$1}}, $leafs[$il] if ($leaf_id =~ /^[']*($fmtGenotype)/); # Sometimes ids have \'
        }
        $debug && print STDERR "\n";
#        $debug && print STDERR "$subname: \$groups=".Dumper($groups)."End of \$groups\n";
        my $group_ids = [ sort {$#{$groups->{$b}} <=> $#{$groups->{$a}}} (keys %$groups) ];
        $debug && print STDERR "$subname: \$focus=".$focus->id." \$group_ids='@$group_ids'\n";

        # Find LCA for each group, these lcas are good choices to re-root the tree
        my $lcas;
        for my $id (@$group_ids) {
          $debug && print STDERR "$subname: \$group_id=$id seqs=$#{$groups->{$id}}\n";
          next if ($#{$groups->{$id}}<1); # Skip those groups with only 1 member
          my $lca = Genotype_newickBI::find_lca( $tree, $groups->{$id});
          push @$lcas, $lca;
          $debug && print STDERR "$subname: for \$group_id=$id \$lca=$lca\n";
        }
#        return if ($#{$lcas}<1); # Quit if there is only 1 group. Can't remember why

        # Use the lca furthest from focus as new root
        my $d = 0.0;
        my $old_root = $tree->get_root_node;
        for my $i (0 ..$#{$lcas}) {
            my $d1 = $tree->distance(-nodes => [ $lcas->[$i], $focus ]);
            $debug && print STDERR "$subname: \$i=$i \$d=$d \$d1=$d1 \$new_root=$new_root\n";
            next if ($d1 <= $d);
            $new_root = $lcas->[$i];
            $d = $d1;
            $debug && print STDERR "$subname: \$i=$i \$d=$d \$d1=$d1 \$new_root=$new_root\n";
        }
#        $debug && print STDERR "$subname: \$new_root=".Dumper($new_root)."End of \$new_root\n";
        $debug && print STDERR "$subname: \$new_root=$new_root \$old_root=$old_root\n";
        $debug && $new_root && print STDERR "$subname: distance=".$tree->distance(-nodes => [ $new_root, $old_root ])."\n";
        if ($new_root && $new_root->ancestor && ($tree->distance(-nodes => [ $new_root, $new_root->ancestor ])>0)) {
            $new_root = $new_root->create_node_on_branch(-fraction => 0.5);
        }
        if ($debug) {
            my $out = new Bio::TreeIO(-file => '>mytree.new', -format => 'newick');
            $out->write_tree($tree);
        }
        $tree->reroot($new_root) if ($new_root && ($new_root != $old_root));
        if ($debug) {
            my $out = new Bio::TreeIO(-file => '>mytree_re_root_1.new', -format => 'newick');
            $out->write_tree($tree);
        }

    }

    $debug && print STDERR "$subname: leaving subroutine\n";
    return;
} # sub re_root_tree


=head2 find_lca

Takes a tree and a number of nodes, Find the LCA for the nodes
Returns undef if no node is provided, or the node if there is only one

=cut

sub find_lca {
    my ($tree, $leafs) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'Genotype_newickBI::find_lca';

    # Find the LCA for the input nodes
    $debug && print STDERR "$subname: \$#leafs=$#{$leafs}\n";
    my $lca;
    return $lca if ($#{$leafs}<0);
    $lca = $leafs->[0];
    return $lca if ($#{$leafs}==0);

    for my $i (1 .. $#{$leafs}) {
          $lca = $tree->get_lca(-nodes => [ $lca, $leafs->[$i] ]);
          ($debug && $i<3) && print STDERR "$subname: for \$i=$i \$lca=$lca\n";
    }

#    $debug && print STDERR "$subname: leaving subroutine\n";
    return $lca;
} # sub find_lca


=head2 find_focalNode

Find the focal node within a phylogeny tree
Takes a tree, and the id of a focal node, returns the focal node

=cut

sub find_focalNode {
    my ($tree, $focalNode) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'Genotype_newickBI::find_focalNode';
    my $dumper_debug = 0;

    $debug && print STDERR "$subname: input \$focalNode=$focalNode \$tree=$tree\n";
    my @node;
    my $ids = [$focalNode, "'$focalNode'"];
    push @$ids, $focalNode =~ s/\'//g;
    $debug && print STDERR "$subname: \@ids=$#{$ids}\n".Dumper($ids)."\n";

    for my $i (0 .. $#{$ids}) {
        $debug && print STDERR "$subname: $i) Trying to find FocalNode=$ids->[$i]\n";
        next if (!$ids->[$i]);
        @node = $tree->find_node(-id => $ids->[$i]);
        next if $#node<0;
        $focalNode = $ids->[$i];
        last;
    }

#    $debug && print STDERR "$subname: \@node=$#node\n".Dumper(@node)."\n";
    &mydie("ERROR: node \$focalNode='$focalNode' not present in \$tree=$tree\n") if $#node==-1;
    &mydie("ERROR: node \$focalNode='$focalNode' not distinct\n") if $#node>0;

    my $focus = $node[0];

    $debug && print STDERR "$subname: found \$focalNode=$focalNode in \$tree=$tree\n";
    $debug && print STDERR "$subname: leaving subroutine\n";
    return $focus;
} # sub find_focalNode


sub mydie() {
    my ($message) = @_;

    my $debug = 0 || $debug_all;
    $message = "Genotype_newickBI::mydie: ".$message;
    print $message . "\n" if $debug;

    `/bin/cp -f $treeFile $failDir` if $PRESERVE_FAILURES;
    `echo $treeFile $focalNode >> $failDir/LOG` if $PRESERVE_FAILURES;
    `echo $message >> $failDir/LOG` if $PRESERVE_FAILURES;

    die $message;
} # sub mydie


1;
