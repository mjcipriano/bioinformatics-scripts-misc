#!/usr/bin/perl

use strict;



my $annotation_file = $ARGV[0];
my $keyword_file = $ARGV[1];


open(KEYS, $keyword_file);

my $keys;

while(<KEYS>)
{
	my $line = $_;
	chomp($line);

	# lower case the whole thing
	$line = lc($line);

	#Convert spaces to dashes
	$line =~ s/\ /-/gi;

	# Remove dashes
	$line =~ s/\-//g;

	$keys->{$line} = 1;
}

my $hits;

open(ANN, $annotation_file);
while(<ANN>)
{
	my $line = $_;
	chomp($line);

	# lower case the whole thing
	$line = lc($line);

	#Convert spaces to dashes
	$line =~ s/\ /-/gi;

	# Remove dashes
	$line =~ s/\-//g;

	# Now search this line with all of the keys

	while(my ($term, undef) = each(%$keys))
	{
		if($line =~ /$term/i)
		{
			# Take the first field which is the uid
			my ($uid) = $line =~ /^(\d+)/;
			$hits->{$uid} = 1;
		}
	}
}

my @results = sort{ $a <=> $b } keys %$hits;

foreach my $hit(@results)
{
	print $hit . "\n";
}


