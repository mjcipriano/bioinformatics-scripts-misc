#!/usr/bin/perl



while(<>)
{
	my $line = $_;
	chomp($line);
	my ($supercontig, $basenumber, $consensusbase, $coverage, $data) = split("\t", $line);
	print join("\t", $supercontig, "PILEUP", "coverage", $basenumber, $basenumber, $coverage, ".", ".", "coverage cdna_coverage ; data \"" . $data . "\" ; ") . "\n";
}
