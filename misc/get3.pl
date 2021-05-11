#!/usr/bin/perl



use Bio::SeqIO;


my $seqio = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'fasta');
my $seqout = Bio::SeqIO->new(-file=>">" . $ARGV[1], -format=>'fasta');

while(my $seq = $seqio->next_seq)
{

	# Get the codon usage of a sliding window
	if($seq->subseq(1,3) ne 'ATG')
	{
		print join("\t", $seq->display_id(), $seq->subseq(1,3)) . "\n";
	} else
	{
		$seqout->write_seq($seq);
	}
}
