#!/usr/bin/perl

use Bio::DB::GenBank;
use strict;


my $gb = new Bio::DB::GenBank;

while(<>)
{
	my $line = $_;
	chomp($line);
	my @accs = split(",", $line);
	foreach my $acc (@accs)
	{
		my $seq = $gb->get_Seq_by_acc($acc); # Accession Number
		foreach my $feat_object ($seq->get_SeqFeatures) 
		{
	        	next unless $feat_object->primary_tag eq "CDS";
			foreach my $tag ($feat_object->get_all_tags) 
			{
				foreach my $value ($feat_object->get_tag_values($tag))
				{
					if($tag eq 'protein_id')
					{
						print $value . "\n";
					}
				}
			}
		}
	}

}


