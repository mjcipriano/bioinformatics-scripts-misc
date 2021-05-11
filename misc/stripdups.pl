#!/usr/bin/perl

use strict;

use Bio::SeqIO;


my %seqhash;
my $in  = Bio::SeqIO->new(-file => 'all.fas' , -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">all_striped.fas" , -format => 'fasta');

while(my $seq = $in->next_seq)
{
	$seqhash{$seq->display_name} = $seq;
}

while(my ($key, $val) = each (%seqhash))
{
	$out->write_seq($val);
}


