#!/usr/bin/perl


use strict;

my $dir = $ARGV[0];

search_dir($ARGV[0]);


sub search_dir
{
	my $subdir = shift;
	opendir(THISDIR, $subdir);
	my @files = readdir(THISDIR);
	close(THISDIR);
	foreach my $file (@files)
	{
		if($file =~ /^\.+$/)
		{
			print $file . "\n";
			next;
		}

		if(-d $dir . "/" . $file)
		{
			search_dir($subdir . "/" . $file);
			# print "$subdir/$file\n";
			next;
		}
		if($file =~ /\.(c|h|cpp)$/)
		{
			#print $file . "\n";
			print_strings_in_file("$subdir/$file");
			
		}
	}

}

sub print_strings_in_file
{
	my $filename = shift;
	open(CFILE, $filename);
	my $file_string = "";
	while(<CFILE>)
	{
		$file_string .= $_;
	}
	# Now search for my regex
	my $match_string;
	while($file_string =~ /(\bui.*?\;)/msg)
	{
		my $start_line = get_line_number($file_string, $-[0]);
		my $end_line = get_line_number($file_string, $+[0]);
		my $match_string = $1;
		while($match_string =~ /(\w+)\((\w+)\)/msg)
		{
			print $filename . "\t" . $start_line . "-" . $end_line . "\t" . $1 . "\t" . $2 . "\n\n";
		}
	}

}

sub get_line_number
{
	my $string = shift;
	my $start = shift;

	my $substart = substr($string, 0, $start);
	my $start_lines = ($substart =~ tr/\n//);
	
	return $start_lines  + 1;
}

