#!/usr/bin/perl



use strict;

my $file = $ARGV[0];
my $type = $ARGV[1];
my $source = $ARGV[2];

open(FILE, $file);
my $id_num = 1;

while(<FILE>)
{
	my $line = $_;
	chomp($line);

	my ($segment, $gsid, $start, $end, $strand, $len, $org, $loc) = split("\t", $line);

	if(!defined($segment))
	{
		next;
	}
	print join("\t", $gsid, $source, $type, $start, $end, "0", $strand, ".", "ID=$type-$id_num") . "\n";
	$id_num++;

}
