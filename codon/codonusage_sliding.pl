#!/usr/bin/perl

use strict;

use Bio::SeqIO;
use Chart::Clicker;
use Chart::Clicker::Data::DataSet;
use Chart::Clicker::Data::Series;

my $window_size = 66;
my $window_movement = 12;
my $start_base = 4;

my $min_num_seqs = 200;
my $graph_width = 1200;
my $graph_height = 600;

my $seqinfile = $ARGV[0];

if(defined($ARGV[1]))
{
	$window_size = $ARGV[1];
}
if(defined($ARGV[2]))
{
	$window_movement = $ARGV[2];
}

my $seqio = Bio::SeqIO->new(-file=>$seqinfile, -format=>'fasta');
my $dir = $window_size . "-" . $window_movement;
mkdir($dir);
mkdir($dir . "/fasta");
mkdir($dir . "/tables");
mkdir($dir . "/summary");
mkdir($dir . "/graphs");

my @seqs;
while(my $seq = $seqio->next_seq)
{
	push(@seqs, $seq);
}

# get overall codon usage
system("cusp $seqinfile -outfile $dir/tables/window_overall.cusp");
system('perl -pi -e "s/^#CdsCount.+\n//s" ' . "$dir/tables/window_overall.cusp");

my $ismore = 1;
my $window_number = 1;
my $this_start = $start_base;
my $this_end = $start_base + $window_size-1;
my $alldata;
my %data;
my %codoninfo;
my $total_seqs = 0;

print "Starting analysis.\n";
while($ismore)
{
	
	$ismore = 0;
	my $seqio = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'fasta');
	my $basename = "window_" . sprintf("%06d",$window_number);
	my $window_fasta_file = $dir . "/fasta/" . $basename . ".fasta";
	my $window_cusp_file = $dir . "/tables/" . $basename . ".cusp";
	my $thiswindow_out = Bio::SeqIO->new(-file=>">". $window_fasta_file, -format=>'fasta');
	my $num_seqs = 0;
	foreach my $seq (@seqs)
	{
		if($seq->length > $this_end)
		{
			$thiswindow_out->write_seq($seq->trunc($this_start, $this_end));
			$num_seqs++;
			$ismore = 1;
		}
	}
	print "$num_seqs sequences in window set $window_number from $this_start-$this_end.\n";;
	if($ismore && $num_seqs > $min_num_seqs )
	{
		system("cusp $window_fasta_file -outfile $window_cusp_file"); 
		system('perl -pi -e "s/^#CdsCount.+\n//s"' . " $window_cusp_file");
		# Slurp in the results
		if(-e $window_cusp_file)
		{
			open(CUSP, $window_cusp_file);
			while(<CUSP>)
			{
				my $line = $_;
				chomp($line);
				if($line =~ /^#/)
				{
					next;
				}
				my ($codon, $aa, $fraction, $freq, $number) = split(/\s+/, $line);
				my %rec;
				$rec{codon} = $codon;
				$rec{aa} = $aa;
				$rec{fraction} = $fraction;
				$rec{freq} = $freq;
				$rec{number} = $number;
				$alldata->{$window_number} = \%rec;
				$codoninfo{$codon} = $aa;
				push(@{$data{fractions}{$codon}}, $fraction);
				push(@{$data{frequency}{$codon}}, $freq);
				push(@{$data{numbers}{$codon}}, $number);

			}
		}
		$window_number++;
		$this_start += $window_movement;
		$this_end += $window_movement;
	} else
	{
		print "$num_seqs below $min_num_seqs minimum, ignoring rest of sequences.\n";
		last;
	}
}

print "Creating Graphs\n";

my @types = ("fractions", "frequency", "numbers");
foreach my $type (@types)
{
	print "Creating $type graph file.\n";
	my $last_aa;
	open(DATAFILE, ">$dir/summary/$type.txt");
	my $graph_format = "png";
	my $cc = Chart::Clicker->new(format=>$graph_format);
	my @ds = ();
	open(HTML, ">$dir/graphs/$type.html");
	print HTML "<html><head><title>$type</title></head><body>";
	foreach my $codon (sort {$codoninfo{$a} cmp $codoninfo{$b} } keys %codoninfo )
	{
		print DATAFILE join("\t", $codoninfo{$codon}, $codon, @{$data{$type}{$codon}}) . "\n";
		
		if(defined($last_aa) && $last_aa ne $codoninfo{$codon})
		{
			if(defined($last_aa))
			{
				my $dataseries =  Chart::Clicker::Data::DataSet->new(series =>\@ds);
				$cc->add_to_datasets($dataseries);
				my $imgfname = $last_aa . "-codons-$type." . $graph_format;
				$cc->write_output($dir . "/graphs/" . $imgfname);
				print HTML "<img src='$imgfname' alt='$last_aa'><br>";
			}
			$cc = Chart::Clicker->new(format=>$graph_format, width=>$graph_width, height=>$graph_height);
			$cc->get_context('default')->domain_axis->format('%d');
			$cc->get_context('default')->domain_axis->label('AA position');
			$cc->get_context('default')->range_axis->label($type);
			if($type eq 'fractions')
			{
				$cc->get_context('default')->range_axis->range->max(1);
				$cc->get_context('default')->range_axis->range->min(0);
				$cc->get_context('default')->range_axis->range->ticks(10);
			}
			# For stacked graph add following
			#$cc->get_context('default')->renderer->additive(1);
			
			@ds = ();
			$cc->title->text("Codon $type for " . $codoninfo{$codon});
		}
		my @size_array = 1 .. scalar @{$data{$type}{$codon}};
		my $codonseries = Chart::Clicker::Data::Series->new(name=>$codon);
		foreach my $ar (@size_array)
		{
			my $pos = int(($window_size/2)/3) + int((($ar-1)*$window_movement)/3);
			#my $pos = (($window_movement * ($ar-1)/3) + 1) + int(($window_size/3)/2);
			#my $pos = int($window_movement * $ar - ($window_size/2)) /3;
			$codonseries->add_pair($pos, $data{$type}{$codon}->[$ar-1]);
		}
		push(@ds, $codonseries);
		$last_aa = $codoninfo{$codon};
	}

	my $dataseries =  Chart::Clicker::Data::DataSet->new(series =>\@ds);
	$cc->add_to_datasets($dataseries);
	my $imgfname = $last_aa . "-codons-$type." . $graph_format;
	$cc->write_output($dir . "/graphs/" . $imgfname);
	print HTML "<img src='$imgfname' alt='$last_aa'><br>";

	close(DATAFILE);
	print HTML "</body></html>";
	close(HTML);

}

