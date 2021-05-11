#!/usr/bin/perl


use strict;

use Bio::DB::Taxonomy;

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');


#my $genomefile = "/data/projects/kinetochore/data2/genomes/genomelist.txt";
my $genomefile = $ARGV[0];

open(GENOMEFILE, $genomefile);

my @taxon_prjids = ();
my %taxon_info;
my @labels = ();
my %label2row;

while(my $line = <GENOMEFILE>)
{
	chomp($line);
	if( $line =~ /^\#/)
	{
		$line =~ s/\#//;
		@labels = split("\t", $line);
		my $lid = 0;
		foreach my $label (@labels)
		{
			$label2row{$label} = $lid;
			$lid++;
		}
	} else
	{
		my @thisarr = split("\t", $line);
		$taxon_info{$thisarr[2]} = \@thisarr;
		push(@taxon_prjids, $thisarr[2]);
	}
}

my %taxoncache = {};

foreach my $prjid (@taxon_prjids)
{
	my $lineprint = '';

	my $taxonid = $taxon_info{$prjid}->[1];
	if(exists($taxoncache{$taxonid}))
	{
		print join("\t", @{$taxon_info{$prjid}}) . "\t";
		print $taxoncache{$taxonid};
		print "\n";

	} else
	{
		my $taxon = $db->get_taxon(-taxonid => $taxonid);
		my @common_names = $taxon->common_names();

		$lineprint .= $taxon->division() . "\t";
		$lineprint .= join(",", @common_names) . "\t";
		
		my $t = $taxon;
		my @all_taxa = ();
		while(defined($t))
		{
			unshift(@all_taxa,  @{$t->name('scientific')});
			my $a = $db->ancestor($t);
			$t = $a;
		}
		$lineprint .= join(",", @all_taxa);
		$taxoncache{$taxonid} = $lineprint;	

		print join("\t", @{$taxon_info{$prjid}}) . "\t";
		print $lineprint;
		print "\n";
	}
}

