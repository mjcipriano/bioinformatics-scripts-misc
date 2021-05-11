#!/usr/bin/perl


while(<>)
{
	my $line = $_;
	$line =~ tr/ATGCatgc/TACGtacg/;
	$line = reverse($line);
	print $line;
}
