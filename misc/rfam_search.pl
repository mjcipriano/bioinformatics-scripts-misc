#!/usr/bin/perl
 
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
 
my $organism = $ARGV[0];

my $rfam_file = $ARGV[1];
my $fasta_input = $ARGV[2];

 
 

my $cmsearch_bin = '/opt/BioBrew/infernal/bin/cmsearch'; 
 
 

# Make a temporary input file

#system("mkdir temp");

#system("$cmsearch_bin $rfam_file $fasta_input > $outfile");



# Open the Rfam directory

my $rfam_dir = "/blastdb/rfam";
my $seed_file = "Rfam.seed";

open(RFAM, $rfam_dir . "/" . $seed_file);

my $rfam_hash;
my $rfam_id_hash;
my $rfam_desc_hash;
my $rfam_wsize_hash;
my $rfam_comment_hash;

my $temp_ac;
my $temp_id;
my $temp_desc;
my $temp_comment;
my $temp_wsize;
my $last_id = '0';

while(<RFAM>)
{
	if($_ =~ /^\#=GF\ AC/)
	{
		($temp_ac) = $_ =~ /^\#=GF\ AC\s+(.+)$/;
		$rfam_hash->{$temp_ac} = 1;
	} elsif($_ =~ /^\#=GF\ ID/)
	{
		($temp_id) = $_ =~ /^\#=GF\ ID\s+(.+)$/;
		$rfam_id_hash->{$temp_ac} = $temp_id;
	} elsif($_ =~ /^\#=GF\ DE/)
	{
		($temp_desc) = $_ =~ /^\#=GF\ DE\s+(.+)$/;
		$rfam_desc_hash->{$temp_ac} = $temp_desc;
	} elsif($_ =~ /^\#=GF\ BM\s+cmsearch/)
	{
		($temp_wsize) = $_ =~ /^\#=GF\ BM\s+cmsearch.+\-W\ (\d+)/;
		$rfam_wsize_hash->{$temp_ac} = $temp_wsize;
	}
	
}
while(my ($key, $val) = each(%$rfam_hash))
{
	print join("\t", $key, $rfam_id_hash->{$key}, $rfam_desc_hash->{$key}, $rfam_wsize_hash->{$key}) . "\n";
}

# Run Search here and redirect to rfam.out

my $outfile = "rfam.out";

my $cm_name = '5S';

open(OUTFILE, $outfile);



my $sequence_name = '';
my $ref_sequence = '';
my $homology_string = '';
my $query_string = '';
my $secondary_string = '';
my $score = '';
my $hit_number = '';
my $query_start = 0;
my $query_end = 0;
my $ref_start = 0;
my $ref_end = 0;
my $line_size = 0;
my $next_sequence_name = '';

my $last_line_type = 'end';
my $first = 0;

while(<OUTFILE>)
{
	my $line = $_;
	if($line =~ /^sequence/)
	{
		($next_sequence_name) = $line =~ /^sequence: (.+)$/;
		
	} elsif($line =~ /^\s$/)
	{
		# This is a whitespace line, do nothing
	} elsif($line =~ /^hit/)
	{
		# We are ready to print out the last hit
		my $direction = "+";
		my $temp_start = $query_start;
		my $temp_end = $query_end;
		if($query_start > $query_end)
		{
			my $temp = $temp_end;
			my $temp_end = $temp_start;
			my $temp_start = $temp;
			$direction = "-";
			
		}
		if($first)
		{
			print_results($cm_name, $query_string, $sequence_name, $query_start, $query_end, $score);
		} else
		{
			$first = 1;
		}
#		print join("\t", $sequence_name, "RFAM", "RNA", $temp_start, $temp_end, $score, $direction, ".", "Sequence name ; " );
		
		print $ref_start . $secondary_string . $ref_end . "\n" . $ref_sequence . "\n" . $homology_string . "\n" . $query_start . $query_string . $query_end . "\n\n";
		$sequence_name = $next_sequence_name;

		
		$hit_number = '';
		$ref_sequence = '';
		$homology_string = '';
		$query_string = '';
		$secondary_string = '';
		$score = '';
		$hit_number = '';
		$query_start = 0;
		$query_end = 0;
		$ref_start = 0;
		$ref_end = 0;
		$line_size = 0;
		$last_line_type = 'hit';

		($hit_number, $query_start, $query_end, $score) = $line =~ /^hit\ (\d+)\s+\:\s*(\d+)\s*(\d+)\s*([0-9\.]+)\ bits\s*$/;

	} elsif($line =~ /^CPU/)
	{		
		print $ref_start . $secondary_string . $ref_end . "\n" . $ref_sequence . "\n" . $homology_string . "\n" . $query_start . $query_string . $query_end . "\n\n";
		if($first)
		{
			print_results($cm_name, $query_string, $sequence_name, $query_start, $query_end, $score);
		} else
		{
			$first = 1;
		}
	} elsif($line =~ /^memory/)
	{
	} else
	{
		if($last_line_type eq 'hit' || $last_line_type eq 'query')
		{
			if($last_line_type eq 'hit')
			{
				$ref_start = 0;
				$ref_end = 0;
			}
			$last_line_type = 'secondary';
			my ($seq_line) = $line =~ /^\s{11}(.+)/;
			$secondary_string .= $seq_line;
			$line_size = length($seq_line);
		} elsif($last_line_type eq 'secondary')
		{
			$last_line_type = 'ref_sequence';
			my ($start, $seq_line, $end) = $line =~ /^\s*(\d+)\s*(.+)\s+(\d+)$/;
			if($ref_start == 0)
			{
				$ref_start = $start;
			}
			if($ref_end == 0)
			{
				$ref_end = $end;
			}
			if($ref_start < $ref_end)
			{
				if($ref_end < $end)
				{
					$ref_end = $end;
				}
			} else
			{
				if($ref_end > $end)
				{
					$ref_end = $end;
				}
			}
			$ref_sequence .= $seq_line;
		} elsif($last_line_type eq 'ref_sequence')
		{
			$last_line_type = 'homology';
			my ($hom_line) = $line =~ /^\s{11}(.+)\s/;
			$homology_string .= $hom_line;
			if(length($hom_line) != $line_size)
			{
				my $num_spaces = $line_size - length($hom_line);
				for(1..$num_spaces)
				{
					$homology_string .= " ";
				}
			}
		} elsif($last_line_type eq 'homology')
		{
			$last_line_type = 'query';
			my ($start, $seq_line, $end) = $line =~ /^\s*(\d+)\s*(.+)\s+(\d+)$/;

			$query_string .= $seq_line;
		}
	}
}

sub print_results
{
	my $cm_name = shift;
	my $sequence = shift;
	my $sequence_name = shift;
	my $start = shift;
	my $stop = shift;
	my $score = shift;

	 $sequence =~ s/\-//gi;

	
	open(FASTA, ">>", $cm_name . '.fas');
	open(GFF, ">>", $cm_name . '.gff');
	print FASTA ">$sequence_name [$start-$stop]\n";
	print FASTA $sequence . "\n";
	my $dir = "+";
	if($stop < $start)
	{
		$dir = "-";
		my $temp = $start;
		$start = $stop;
		$stop = $temp;
	}
	print GFF join("\t", $sequence_name, 'rfam', 'misc_RNA', $start, $stop, $score, $dir, '.', 'rfam ' . $cm_name) . "\n";

}




sub create_alignment
{
	my $filename = shift;
	my $fasta_file = shift;

}
