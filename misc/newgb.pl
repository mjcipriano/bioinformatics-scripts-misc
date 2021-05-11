#!/usr/bin/perl



use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use strict;

my $in = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'genbank');
my $out = Bio::SeqIO->new(-file=>">" . $ARGV[1], -format=>'genbank');

my $seq_add = "acgtacatcgtttcaaatgaagtactcaatgtgcactattctgggggatacattataatctttactatggagggga
agatatggcggttaaaccttctcgatgacaagtcacttgaaaactaaaaaatcactgagcagccggatggatgggctcac
acagccgttcaccgaccagccaatggaggtcatgaccgatgagcagctcaaggaaatacatgaaaagaagatcatcgacg
caatcaacggctccatttttggacttgccaaagtggtcaaggaaactcagcacatcgaagtcccagtggccactctggat
caggggttgcctctgaagccgtcaagcgtccttgcaacccttatgatgttctacaacacagcatatctccttcccttgtc
caaagttatcgatgacatccccccgcaatttaaagcgcttttcagcgagcattacgccgaaaagtatcaatatatcctca
atatctatgctcaatacggtattttgaaagagttcctagccgattgcaagtatccgggcttttcctttactgacatcttt
atgcctgcattcaagcgctacaaagcccagctgcacaatatctctgaatttttatggttccgtgaggcggtaattgaaat
ccttcactcagccctaacagaatacaaaactgcgtacactcaaaataagcctgttttggaatctgcatcagcccttaagc
taggaattagtcgcataacagatgaaactgtaaagctccgggagacggcctctcagctgccatcgcggatcgatcccctt
attcaacgatacaacagagcactcgagagccaacaacgtcttgataaggaaaatagcgagttcaccaaccggatacttgt
gctgactaacgaactagttgagctggaaggtcaagtcaagaagaaaaccgaagaatatgaatctgaaaaaagcactctaa
tcactaaaaaacaaaacgatgccgttcttcaggacatacaagctgctgaggccttaggaaatgacctgagaactcagtgc
gatgctctcttatcatcgatcaatatttttgagcagcacaacaatgatcttctggggataaccgacaagataaccgaact
gattgacacagctgatcaatatagtagctccagtgacaatctgcaggcacagttggagagcctacaggtaatgaaggata
ccgttaataccctggacgcgcagattaatgcatccaaacaagaaatagagaatatgaacaatgcaatgcagtatggatct
gccgtagagcaagaacaggacaaccctagtcttcgcgcagatgaagaggagcttcggcgactcaccgaagattgcgaagg
cctggaagccacatatgctaagaagctcagtcaggttagagccctacagcaagaaatcaatgacgagaatgcaaaaataa
ccagggtgagagtggataccgctgaaaagatcgaccagctgtctgacaaggccaatgaaatgatgtcaatgtttgtagcg
caactcacaggaatgataagcaagctcagctccatg";
$seq_add =~ s/\n//g;

my $seq = $in->next_seq();

my $seq_end = $seq->length();

$seq->seq($seq->seq() . $seq_add);
my $new_seqend = $seq->length();
my @feats = ();
my $fiveclipsize = 100;
push(@feats, get_feat($seq_end+1, $new_seqend, 0, 'misc_feature', 'vnti', 'TOPO insert', 'TOPO insert', 21) );
push(@feats, get_feat($seq_end+1,  $seq_end + $fiveclipsize, 0, "5'clip", 'vnti', "gene 5' Clip", "gene 5' Clip", 51) );
push(@feats, get_feat($seq_end + $fiveclipsize + 1, $new_seqend, 0, "CDS", 'vnti', "ORFid ", "Orfid ", 4) );

$seq->add_SeqFeature(@feats);

$out->write_seq($seq);





sub get_feat
{
	my ($start, $end, $strand, $primary, $source_tag, $display_name, $label, $vntifkey) = @_;

	my $seqfeat = Bio::SeqFeature::Generic->new(
						-start=>	$start,
						-end=>		$end,
						-strand=>	$strand,
						-primary=>	$primary,
						-source_tag=>	$source_tag,
						-display_name=>	$display_name,
						-tag=>		{ 	label=>$label,
									vntifkey=>$vntifkey
								}
		);
	return $seqfeat;

}

