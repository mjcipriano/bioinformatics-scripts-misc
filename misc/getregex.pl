#!/usr/bin/perl


use strict;


use Bio::SeqIO;


my $file = $ARGV[0];
my $type = $ARGV[1];
my $source = $ARGV[2];

my $single_dir = 0;

my $seqs = Bio::SeqIO->new(-file=>$file, -format=>"fasta");

my $num = 1;

while(my $seq = $seqs->next_seq)
{
	my $name = $seq->display_id();
	my $seqf = $seq->seq();
	while( $seqf =~ /(GG.{19}.GG)/gi)
	{
		print join("\t", $name, $source, $type, @-[0]+1, @+[0]+1, "0", "+", ".", "ID=$type-$num;Note=$1" ) . "\n";
		$num++;
	}

	if($single_dir)
	{
		my $seqr = $seq->revcom()->seq();

		my $length = $seq->length();

		while( $seqr =~ /(GG.{19}.GG)/gi)
		{
			print join("\t", $name, $source, $type, $length-@-[0], $length+@-[0], "0", "-", ".", "ID=$type-$num;Note=$1" ) . "\n";
			$num++;
		}
	}
}
