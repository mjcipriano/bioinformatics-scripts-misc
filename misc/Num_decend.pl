#!/usr/bin/perl

use Bio::TreeIO;
use strict;



my $treeio = Bio::TreeIO->new(-format => 'nexus' , -file => 'Rhodopsinmeta2.tre');
my $treeout = Bio::TreeIO->new(-format => 'nexus' , -file => '>Metno.tre');
my $treeout2 = Bio::TreeIO->new(-format => 'nexus' , -file => '>orig.tre');

my $tree = $treeio->next_tree;
#my $rootnode = $tree->get_root_node;

$treeout2->write_tree($tree);

foreach my $node ( $tree->get_nodes() ) 
{
	
	if( !defined($node->ancestor))
	{
		print "SKIPPING " . $node->id() . "\n";
		next;
	}
	my $count = count_all_leaves($node);
	print $node->id . "\t" . $node->internal_id . "\t" . $count . "\n";
}

$treeout->write_tree($tree);

sub count_all_leaves
{
	my $node = shift;
	if ($node->is_Leaf)
	{
		$node->id($node->id());
		$node->description($node->id());
		return 0;
	}
	my $count = 0;
	my $keep = 0;
	my $name = $node->id();
	foreach my $innode ( $node->get_all_Descendents) 
	{
		if ( $innode->is_Leaf )
		{
			$count += 1;
			if ( $innode->id =~ m/\>JCVIPEP/   )
			{
			} else
			{
				#print $innode->id() . "\n";
				$keep = 1;
			}
		}

	}
	if ($keep )
	{
		$node->id($name);	
	} else
	{
		
		 $node->remove_all_Descendents();
		my $newname = $node->internal_id() . "_count_" . $count;
		$node->id($newname);
		$node->description($newname);
		#my $newnode = Bio::Tree::Node->new(-id=>$node->internal_id() . "_count_" . $count, -description=>'collapsed_node');
		#$node->add_Descendent($newnode);
		
	}
	
	return $count;

}
