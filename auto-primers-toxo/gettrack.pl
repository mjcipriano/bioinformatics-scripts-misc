#!/usr/bin/perl



use strict;
my $filename = $ARGV[0];
open(FILE, $filename);
# http://toxodb.org/cgi-bin/gbrowse/toxodb/?plugin=TrackDumper;plugin_action=Go;name=TGME49_chrV:60000..70000;view_start=60000;view_stop=319999;plugin_button=Cancel&plugin_button=Configure&plugin_button=Go&TrackDumper.mode=selected&TrackDumper.version=3&TrackDumper.region=All&TrackDumper.disposition=view&plugin_button=Cancel&plugin_button=Configure&plugin_button=Go

my $region_size = 200000;

while(<FILE>)
{
	my $line = $_;
	chomp($line);
	my ($name, $source, $type, $start, $end, $score, $dir, $frame, @rest) = split("\t", $line);
	if($type ne "supercontig" ) # TODO add chromosome
	{	
		next;
	}	
	print $type . " $name\n";
	print $start . "\n";
	my $t_start = $start;
	my $t_end = $start + $region_size -1;
	if($t_end > $end)
	{
		$t_end = $end;
	}

	while($t_start <= $end)
	{
		my $link = join("", 
					"http://toxodb.org/cgi-bin/gbrowse/toxodb/?plugin=TrackDumper;plugin_action=Go;plugin_button=Cancel&plugin_button=Configure&plugin_button=Go&TrackDumper.mode=selected&TrackDumper.version=3&TrackDumper.region=All&TrackDumper.disposition=view&plugin_button=Cancel&plugin_button=Configure&plugin_button=Go;",
					"name=$name:$t_start..$t_end;",
					"view_start=$t_start;view_stop=$t_end;"
				);
		print $link . "\n";
		$t_start = $t_end + 1;
		$t_end = $t_end + $region_size;
		if($t_end > $end)
		{
			$t_end = $end;
		}
	}
	
}
