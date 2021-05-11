#!/usr/bin/perl


use strict;

my $numrows = 119;
my @rows;
my $row = 0;
my $column = 0;
my %values;

my $name;
my $lastname;
my $well;
my $typenum;
my $val;
my $mean;
my $lastmean;

while(<>)
{
	my $line = $_;

	if($line =~ /^\s*$/ || $line =~ /^Group:/ || $line =~ /^Sample/ || $line =~ /\~/ || $line =~ /(Smallest|Largest|R - |SD = |Mean Adj)/)
	{
		next;
	}
	chomp $line;
	$lastname = $name;
	($name, $well, $typenum, $val, $mean) = split("\t", $line);
	if($name eq "")
	{
		$name = $lastname;
	}
	if($mean eq "")
	{
		$mean = $lastmean;
	}
	$values{$name}{$well}{$column} = $val;
	# $rows[$row][$column] = $line;
	$rows[$row][$column] = join("\t", $name, $well, $typenum, $val, $mean);
	
	# print " $row $column is $line \n";
	if($row < $numrows)
	{
		$row++;
	} else
	{
		$row = 0;
		$column++;
	}

}
my $firstrow = 1;
for my $i (0..$column)
{
	if($firstrow)
	{
		print join("\t", $i, $i, $i);
		$firstrow = 0;
	} else
	{
		print "\t" . $i; 
	}
}
print "\n";

my %summary;
foreach my $thisrow (@rows)
{
	my $col = 0;
	$firstrow = 1;
	foreach my $thiscol(@{$thisrow})
	{
		my ($n, $w, $t, $v, $m) = split("\t", $thiscol);
		if($firstrow)
		{
			print join("\t", $n, $w, $v);
			$firstrow = 0;
		} else
		{
			print "\t" . $v;
		}
		$col++;
	}
	print "\n";
}

