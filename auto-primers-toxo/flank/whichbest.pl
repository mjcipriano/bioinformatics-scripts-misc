#!/usr/bin/perl


use strict;

my %allenzymes;

while(<>)
{
	my $line = $_;
	chomp($line);
	my ($have, $dont) = split('\|', $line);
	if($dont =~ /^$/)
	{
		$allenzymes{'none'}++;
		next;
	}
	my @enzymes = split(",", $dont);
	foreach my $enz (@enzymes)
	{
		my ($name, $pos) = split(":", $enz);
		if($name =~ /^$/)
		{
			next;
		}
		$allenzymes{$name}++;	
		# print join(":", $name, $pos) . "\n";
	}
}


while(my ($name, $val) = each(%allenzymes))
{
	print join(" : ", $val, $name) . "\n";
}
