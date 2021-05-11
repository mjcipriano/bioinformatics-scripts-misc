#!/usr/bin/perl
use Bio::SeqIO; 

$in = Bio::SeqIO->new(-file=>'ts_nt.fasta', -format=>'fasta');
$out=Bio::SeqIO->new(-file=>'>ts_nt_n.fasta', -format=>'fasta');

while(my $seq = $in->next_seq)
{
	my ($orfid) = $seq->display_id =~ /^(\d+)/;
	$seq->display_id($orfid);
	$seq->description("");
	$out->write_seq($seq);
}
