#!/usr/bin/perl


use strict;

use Bio::SeqIO;

foreach my $file (@ARGV)
{
	my $gb = $file;


	my $seqio = Bio::SeqIO->new(-format=>'genbank', -file=>$gb);


	while(my $seq = $seqio->next_seq)
	{
		my $lt = '';
		my $cds = '';

		# add Accession isolate country PUBMED JOURNAL
		my $isolate = '';
		my $country = '';
		my $pubmed = '';
		my $journal = '';
		my @refs = $seq->get_Annotations('reference');
		#my $pubmed = join("|", map{$_->pubmed}@refs);
		foreach my $ref (@refs)
		{
			if($ref->pubmed)
			{
				$pubmed = $ref->pubmed();
				($journal) = $ref->location() =~ /(.+),/;
			}
		}
		#my $journal = join("|", map{$_->location() =~ /(.+),/ }@refs);
		#my $journal = join("|", map{$_->location() }@refs);

		foreach my $feature ($seq->get_SeqFeatures)
		{
			if($feature->primary_tag eq 'source')
			{
				$isolate = $feature->has_tag('isolate') ? join("", $feature->get_tag_values('isolate')) : "";
				$country = $feature->has_tag('country') ? join("", $feature->get_tag_values('country')) : "";
			}
			if($feature->primary_tag eq 'CDS')
			{
				print join("\t", 
						$seq->accession(),
						$isolate,
						$country,
						$pubmed,
						$journal,
						# defined($seq->species()) ? $seq->species()->binomial() : "",
						$feature->primary_tag, $feature->start(), $feature->end(), $feature->strand(), 
						$feature->has_tag('locus_tag') ? $feature->get_tag_values('locus_tag') : "", 
						$feature->has_tag('product') ? $feature->get_tag_values('product') : "", 
						$feature->has_tag('source') ? $feature->get_tag_values('source') : "",
				) . "\n";
			}
		}	
	}

}


