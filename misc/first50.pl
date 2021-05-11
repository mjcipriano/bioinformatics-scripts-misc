#!/usr/bin/perl

use strict;

use Bio::SeqIO;



my $in = $ARGV[0];
my $out = $ARGV[1];
my $seqio = Bio::SeqIO->new(-file=>$in, -format=>'fasta');
my $seqout = Bio::SeqIO->new(-file=>">$out", -format=>'fasta');

while(my $seq = $seqio->next_seq)
{
	my $subseq = $seq->trunc(1,45);
	$seqout->write_seq($subseq);
	
}
