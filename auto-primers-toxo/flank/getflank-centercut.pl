#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Restriction::Analysis;

use strict;

my $fasta_file = $ARGV[0];
my $gff_file = $ARGV[1];

open(GFF, $gff_file);

my %chr2seq;

my $seqs = Bio::SeqIO->new(-file=>$fasta_file, -format=>'fasta');
my @gb_slic_files = ( "p3xha-lic-hxg-insert.gb", "pcfp-lic-dhfr-insert.gb", "pmcherry-lic-dhfr-insert.gb", "ptap-lic-hxg-insert.gb", "ptdtomato-lic-hxg-insert.gb", "pyfp-lic-dhfr-insert.gb", "p3xmyc-lic-hxg-insert.gb", "pcfp-lic-hxg-insert.gb", "pmcherry-lic-hxg-insert.gb", "ptdtomato-lic-dhfr-insert.gb", "pyfp-lic-cat-insert.gb", "pyfp-lic-hxg-insert.gb" );

my $slic_files_dir = "slic";
my $slic_out_dir = "/home/mcipriano/gb";
#my $slic_out_dir = "created";

while(my $seq = $seqs->next_seq)
{
	$chr2seq{$seq->display_id} = $seq;
	#print $seq->display_id . "\n";
}

my %cds2coord;
my %gene2dir;
my %gene2chr;

my $default_num_flank = 100;
my $re_digest_min_flank = 250;
my $max_gene_size_slic = 5000; # NOTE: this must be one of the static sizes as well
my @static_sizes = (500,600,700,800,900,1000,1200,1500,2000,3000,4000,5000,6000,8000,10000);

my $target_tm = 65;
my $fosmid_homology_size = 50;
my $re_have_file = "re_have.txt";
my $re_avoid_file = "re_avoid.txt";
my @re_avoid = ( );
my @re_have = ();
my %re_have;
my %re_avoid;


if(defined($re_have_file) && -e $re_have_file)
{
	open(RE, $re_have_file);
	while(<RE>)
	{
		my $line = $_;
		chomp($line);
		push(@re_have, $line);
	}
}

if(defined($re_avoid_file) && -e $re_avoid_file)
{
	open(RE, $re_avoid_file);
	while(<RE>)
	{
		my $line = $_;
		chomp($line);
		push(@re_avoid, $line);
	}
}


map {$re_have{$_} = 1} @re_have;
map {$re_avoid{$_} = 1} @re_avoid;


print join("\t", "gene", "name", "start", "end", "dir", "flank_upstream_seq", "gene_start_seq", "gene_end_seq", "flank_downstream_seq", "fullgene_f", "Tm_fullgene_f", "fullgene_r", "Tm_fullgene_r", "max_cutter_name_pos_size", "cut_locs", "LIC_f", "LIC_tm_f", "LIC_r_nostop", "LIC_tm_r_nostop", "slic_seq_amp_length", "slic_seq_amp_seq", "fosmid_5prime_f", "fosmid_5prime_r", "fosmid_3prime_f", "fosmid_3prime_r", "warnings") . "\n";

while(<GFF>)
{
	my $line = $_;
	chomp($line);
	my ($name, $source, $type, $start, $end, $score, $dir, $frame, @rest) = split("\t", $line);

	my ($gene) = @rest[0] =~ /(TGME49_\d+)/;
	$gene2chr{$gene} = $name;
	if($type eq "CDS")
	{
		push(@{$cds2coord{$gene}}, $start, $end);
		$gene2dir{$gene} = $dir;
	}

		
	
}
# gene, min, max, flank_upstream, gene_start, gene_end, flank_downstream
while(my ($gene, $dir) = each(%gene2dir))
{
	my $name = $gene2chr{$gene};
	my $seq;
	my $re_seq;
	my @sorted_coords = sort {$a <=> $b} @{$cds2coord{$gene}};
	my $min = $sorted_coords[0];
	my $max = $sorted_coords[-1];
	my $start = $min;
	my $end = $max;
	my $flank_upstream = '';
	my $gene_start = '';
	my $gene_end = '';
	my $flank_downstream = '';
	my $outofbounds = 0;
	my $warnings = '';
	my $num_flank = $default_num_flank;

	if($start - $num_flank < 1 )
	{
		$num_flank = $start;
		$warnings .= "short 5prime flank of $num_flank bp size;";
		$outofbounds = 1;
	}
	if($end+$num_flank-1 > $chr2seq{$name}->length())
	{
		if($num_flank < $chr2seq{$name}->length() - $end)
		{
		} else
		{
			$num_flank = $chr2seq{$name}->length() - $end;
			$warnings .= "short 3prime flank of $num_flank bp size";
		}
		$outofbounds = 1;
	}

	if(!exists($chr2seq{$name}))
	{
		# print join("\t", $gene, $name, $start, $end, "nosource") . "\n";
		next;
	}

	if($outofbounds)
	{
		print join("\t", $gene, $name, $start, $end, $dir, $warnings) . "\n";
		next;
	}


	if($dir eq "-")
	{
		$flank_upstream = $chr2seq{$name}->trunc($end+1, $end+$num_flank-1)->revcom();
		$gene_start = $chr2seq{$name}->trunc($end-$num_flank, $end)->revcom();
		$gene_end = $chr2seq{$name}->trunc($start, $start+$num_flank-1)->revcom();
		$flank_downstream = $chr2seq{$name}->trunc($start-$num_flank, $start-1)->revcom();
		$seq = $chr2seq{$name}->trunc($start, $end)->revcom();

	} elsif($dir eq "+")
	{
		$flank_upstream = $chr2seq{$name}->trunc($start-$num_flank, $start-1);
		$gene_start = $chr2seq{$name}->trunc($start, $start+$num_flank-1);
		$gene_end = $chr2seq{$name}->trunc($end-$num_flank, $end);
		$flank_downstream = $chr2seq{$name}->trunc($end+1, $end+$num_flank-1);
		$chr2seq{$name}->trunc($start, $end);
		$seq = $chr2seq{$name}->trunc($start, $end);

	} else
	{
		next;
	}

	if($seq->seq() =~ /N/)
	{
		$warnings .= "N's within gene sequence;";
	}
	if($flank_downstream->seq() =~ /N/)
	{
		$warnings .= "N's within downstream sequence;";
	}

	if($flank_upstream->seq() =~ /N/)
	{
		$warnings .= "N's within upstream sequence;";
	}
	# Now find out primers which are a specific tm
	my ($fullgene_f, $fullgene_r) = get_primers($seq->seq(), 20, $target_tm);




	############################## START SLIC ###################################################

	my $max_cutter_name = "none";
	my $max_cutter_slic_pos;
	my $fragment_size;
	my $cut_locs_unq;
	my $cut_locs_unq_nothave;
	my $slicgene_f;
	my $slicgene_r;
	my $slicgene_r_nostop;
	my $slic_seq_amp;
	my $slic_seq_amp_nostop;
	my $previous_warnings = $warnings;
	my $slic_warnings = '';
	

	foreach my $gene_size_slic (@static_sizes)
	{
		if($gene_size_slic > $max_gene_size_slic)
		{
			last;
		}

		if($dir eq "-")
		{
			if($start + $gene_size_slic > $chr2seq{$name}->length())
			{
				$re_seq = $chr2seq{$name}->trunc($start, $chr2seq{$name}->length())->revcom();
			} else
			{
				$re_seq = $chr2seq{$name}->trunc($start, $start + $gene_size_slic)->revcom();
			}

		} elsif($dir eq "+")
		{

			if($end - $gene_size_slic < 1)
			{
				$re_seq = $chr2seq{$name}->trunc(1, $end);
			} else
			{
				$re_seq = $chr2seq{$name}->trunc($end - $gene_size_slic, $end);
			}

		}





		# Now find where the best unique restriction site is on the cterm for single cut slic cloning
		#my $slic_seq = $seq;
		#if($seq->length > $gene_size_slic)
		#{
		#	$slic_seq = $seq->trunc($slic_seq->length() - $gene_size_slic+1, $slic_seq->length)
		#}
		$slic_warnings = '';
		my $slic_seq = $re_seq;
		my $ra = Bio::Restriction::Analysis->new(-seq=>$slic_seq);
		my $single_cutters = $ra->cutters(1);
		my $cut_locs = '';
		my @cut_locs;
		my $max_cutter_pos = 1;
		$max_cutter_name = 'none';
		my $gene_region_size = $slic_seq->length();
		map { 
			if(exists($re_have{$_->name}) && !exists($re_avoid{$_->name})) 
			{
				$cut_locs .= "," . join(":", $_->name, join(",", $ra->positions($_->name)));
				# This only works for single cutters
				my $this_cutter_pos = join("", $ra->positions($_->name));
				if( 	($this_cutter_pos > $max_cutter_pos) && 
					($gene_region_size - $this_cutter_pos >= $re_digest_min_flank) && 
					($this_cutter_pos >= $re_digest_min_flank)
				)
				{
					$max_cutter_pos = $this_cutter_pos;
					$max_cutter_name = $_->name;
				}
			}
		    } $single_cutters->each_enzyme;

		my $slic_pos = $max_cutter_pos - $re_digest_min_flank+1;
		if($slic_pos < 1)
		{
			$slic_pos = 1;
		}

		my $slic_closest_size = closest_size(\@static_sizes, $slic_seq->length() - $slic_pos);
		if($slic_closest_size > $gene_size_slic || $slic_closest_size < 0)
		{
			# Do nothing, this will be an odd size
		} else
		{
			$slic_pos = $slic_closest_size;
		}
		# Get the position of this portion of the gene
		$slic_seq_amp = undef;
		$slic_seq_amp_nostop = undef;
		# TODO something better than this

		if($slic_pos < 1 )
		{
			$slic_warnings .= "no slic within parameters, amplify gene and hope to find something within this;";
			$slic_pos = $end - $start+1;
		}

		if($dir eq "-")
		{
			if($start + $slic_pos >  $chr2seq{$name}->length() )
			{
				$slic_pos = $chr2seq{$name}->length() - $start;
				$slic_warnings .= "slic flank position past end of contig, using end of contig for pcr, check that re digest is within this region;";
			}
			$slic_seq_amp = $chr2seq{$name}->trunc($start, $start+$slic_pos-1)->revcom();

		} elsif($dir eq "+")
		{
			if($end-$slic_pos < 1)
			{
				$slic_pos = $end;
			}
			$slic_seq_amp = $chr2seq{$name}->trunc($end-$slic_pos+1, $end);
		}
		$slic_seq_amp_nostop = $slic_seq_amp->trunc(1, $slic_seq_amp->length()-3);
		($slicgene_f, $slicgene_r) = get_primers($slic_seq_amp->seq(), 20, $target_tm);
		($slicgene_f, $slicgene_r_nostop) = get_primers($slic_seq_amp_nostop->seq(), 20, $target_tm);

		# Now get new all unique cutters
		
		my $ra_unq = Bio::Restriction::Analysis->new(-seq=>$slic_seq_amp_nostop);
		my $single_cutters_unq = $ra_unq->cutters(1);
		$max_cutter_slic_pos = join(" ", $ra_unq->positions($max_cutter_name));

		$fragment_size = $slic_seq_amp_nostop->length() - $max_cutter_slic_pos +1;
		$cut_locs_unq = '';
		$cut_locs_unq_nothave = '';
		map { 
			if(exists($re_have{$_->name}) && !exists($re_avoid{$_->name})) 
			{
				$cut_locs_unq .= "," . join(":", $_->name, join(",", $ra_unq->positions($_->name)));
				# This only works for single cutters
			} elsif(!exists($re_avoid{$_->name}))
			{
				$cut_locs_unq_nothave .= "," . join(":", $_->name, join(",", $ra_unq->positions($_->name)));
			}
		    } $single_cutters_unq->each_enzyme;

		if($max_cutter_name ne 'none')
		{
			last;
		}

	}
	$warnings = $previous_warnings . $slic_warnings;
	################################ END SLIC #############################################################

	my $fosmid_5prime_f = $flank_upstream->trunc($flank_upstream->length()-$fosmid_homology_size, $flank_upstream->length());
	my $fosmid_5prime_r = $gene_start->trunc(1, $fosmid_homology_size)->revcom();	
	my $fosmid_3prime_f = $gene_end->trunc($gene_end->length()-$fosmid_homology_size, $flank_upstream->length()); # If want to get rid of stop, remove last 3 from this one
	my $fosmid_3prime_r = $flank_downstream->trunc(1, $fosmid_homology_size)->revcom();


	print join("\t", $gene, $name, $start, $end, $dir, $flank_upstream->seq, $gene_start->seq, $gene_end->seq, $flank_downstream->seq, $fullgene_f, int(0.5+Tm($fullgene_f)), $fullgene_r, int(0.5+Tm($fullgene_r)), "$max_cutter_name:$max_cutter_slic_pos:$fragment_size", $cut_locs_unq . "|" . $cut_locs_unq_nothave, $slicgene_f, int(0.5+Tm($slicgene_f)), $slicgene_r_nostop, int(0.5+Tm($slicgene_r_nostop)), $slic_seq_amp_nostop->length(), undef, $fosmid_5prime_f->seq, $fosmid_5prime_r->seq, $fosmid_3prime_f->seq, $fosmid_3prime_r->seq, $warnings) . "\n";


	# TODO remove this and create files
	# next;

	# Create genbank files
	foreach my $gb_file (@gb_slic_files)
	{
		my $gb_seqs = Bio::SeqIO->new(-file=>"$slic_files_dir/$gb_file", -format=>'genbank');

		my $gb_out = Bio::SeqIO->new(-file=>">$slic_out_dir/$gb_file-$gene.gb", -format=>'genbank');
		while(my $gb_seq = $gb_seqs->next_seq)
		{
			my $seq_end = $gb_seq->length();
			$gb_seq->seq($gb_seq->seq() . $slic_seq_amp_nostop->seq());
			my $new_seqend = $gb_seq->length();
			my @feats = ();

			push(@feats, get_feat($seq_end+1, $new_seqend, 0, 'misc_feature', 'slic', 'insert', "slic insert $gene cut with $max_cutter_name", 21) );
			push(@feats, get_feat($seq_end+1, $seq_end+1+length($slicgene_f), 0, 'primer', 'slic', 'primer', "slic_primer_f_tm_" . int(0.5+Tm($slicgene_f)), 21) );
			push(@feats, get_feat($new_seqend - length($slicgene_r_nostop), $new_seqend, 0, 'primer', 'slic', 'primer', "slic_primer_r_tm_" . int(0.5+Tm($slicgene_r_nostop)), 21) );
			$gb_seq->add_SeqFeature(@feats);
			$gb_seq->display_id($gb_seq->display_id . "-$gene");
			$gb_out->write_seq($gb_seq);

		}
	}

}

# TODO prep fosmid clones and print out genbank files for all of these


sub closest_size
{
	my $array = shift;
	my $num = shift;

	foreach my $tnum (sort {$a <=> $b} @{$array})
	{
		if($num <= $tnum)
		{
			return $tnum;
		}
	}
	return -1;

}









=head2 get_primers

=cut

sub get_primers
{
	my $sequence = shift;
	my $size = shift;
	my $by_tm = shift;

	my $f_primer = undef;
	my $r_primer = undef;
	

	# If only by size, return forward and reverse substrings

	if(!defined($by_tm) || $by_tm == 0)
	{
		$f_primer = substr($sequence, 0, $size);
		$r_primer = substr($sequence, -$size);
		$r_primer = revcomp($r_primer);
		
		return ($f_primer, $r_primer);
	}
	
	# Get the forward primer
	my $cur_tm = 0;
	my $cur_size = $size - 1;
	my $last_tm = 0;
	if(!defined($cur_size))
	{
		$cur_size = 0;
	}
	while($cur_tm < $by_tm && $cur_size < 30)
	{
		$cur_size++;
		$last_tm = $cur_tm;
		$cur_tm = Tm(substr($sequence, 0, $cur_size));
	}

	if( abs($last_tm - $by_tm) < abs($cur_tm - $by_tm) )
	{
		$cur_size--;
	}
	$f_primer = substr($sequence, 0, $cur_size);


	# Get the reverse primer by tm
	$cur_tm = 0;
	$cur_size = $size - 1;
	$last_tm = 0;
	if(!defined($cur_size))
	{
		$cur_size = 0;
	}
	while($cur_tm < $by_tm && $cur_size < 30)
	{
		$cur_size++;
		$last_tm = $cur_tm;
		$cur_tm = Tm(revcomp(substr($sequence, -$cur_size)));
	}

	if( abs($last_tm - $by_tm) < abs($cur_tm - $by_tm) )
	{
		$cur_size--;
	}

	$r_primer = substr($sequence, -$cur_size);
	$r_primer = revcomp($r_primer); # reverse complement the primer

	return ($f_primer, $r_primer);

}



=head2 Tm()

 Title   : Tm()
 Usage   : $tm = $primer->Tm(-salt=>'0.05')
 Function: Calculates and returns the Tm (melting temperature) of the primer
 Returns : A scalar containing the Tm.
 Args    : -salt set the Na+ concentration on which to base the calculation.
           (A parameter should be added to allow the oligo concentration to be set.)
 Notes   : Calculation of Tm as per Allawi et. al Biochemistry 1997 36:10581-10594.  Also see
           documentation at http://biotools.idtdna.com/analyzer/ as they use this formula and 
           have a couple nice help pages.  These Tm values will be about are about 0.5-3 degrees 
           off from those of the idtdna web tool.  I don't know why.
=cut

sub Tm  {
    my $sequence = shift;
    
    my $salt_conc = 0.05; #salt concentration (molar units)
    my $oligo_conc = 0.00000025; #oligo concentration (molar units)
 #   if ($args{'-salt'}) {$salt_conc = $args{'-salt'}} #accept object defined salt concentration
    #if ($args{'-oligo'}) {$oligo_conc = $args{'-oligo'}} #accept object defined oligo concentration
    my $length = length($sequence);
    my @dinucleotides;
    my $enthalpy;
    my $entropy;
    #Break sequence string into an array of all possible dinucleotides
    while ($sequence =~ /(.)(?=(.))/g) {
        push @dinucleotides, $1.$2;
    }
    #Build a hash with the thermodynamic values
    my %thermo_values = ('AA' => {'enthalpy' => -7.9,
                                  'entropy'  => -22.2},
                         'AC' => {'enthalpy' => -8.4,
                                  'entropy'  => -22.4},
                         'AG' => {'enthalpy' => -7.8,
                                  'entropy'  => -21},
                         'AT' => {'enthalpy' => -7.2,
                                  'entropy'  => -20.4},
                         'CA' => {'enthalpy' => -8.5,
                                  'entropy'  => -22.7},
                         'CC' => {'enthalpy' => -8,
                                  'entropy'  => -19.9},
                         'CG' => {'enthalpy' => -10.6,
                                  'entropy'  => -27.2},
                         'CT' => {'enthalpy' => -7.8,
                                  'entropy'  => -21},
                         'GA' => {'enthalpy' => -8.2,
                                  'entropy'  => -22.2},
                         'GC' => {'enthalpy' => -9.8,
                                  'entropy'  => -24.4},
                         'GG' => {'enthalpy' => -8,
                                  'entropy'  => -19.9},
                         'GT' => {'enthalpy' => -8.4,
                                  'entropy'  => -22.4},
                         'TA' => {'enthalpy' => -7.2,
                                  'entropy'  => -21.3},
                         'TC' => {'enthalpy' => -8.2,
                                  'entropy'  => -22.2},
                         'TG' => {'enthalpy' => -8.5,
                                  'entropy'  => -22.7},
                         'TT' => {'enthalpy' => -7.9,
                                  'entropy'  => -22.2},
                         'A' =>  {'enthalpy' => 2.3,
                                  'entropy'  => 4.1},
                         'C' =>  {'enthalpy' => 0.1,
                                  'entropy'  => -2.8},
                         'G' =>  {'enthalpy' => 0.1,
                                  'entropy'  => -2.8},
                         'T' =>  {'enthalpy' => 2.3,
                                  'entropy'  => 4.1}
                        );
    #Loop through dinucleotides and calculate cumulative enthalpy and entropy values
    for (@dinucleotides) {
       $enthalpy += $thermo_values{$_}{enthalpy};
       $entropy += $thermo_values{$_}{entropy};
    }
    #Account for initiation parameters
    $enthalpy += $thermo_values{substr($sequence, 0, 1)}{enthalpy};
    $entropy += $thermo_values{substr($sequence, 0, 1)}{entropy};
    $enthalpy += $thermo_values{substr($sequence, -1, 1)}{enthalpy};
    $entropy += $thermo_values{substr($sequence, -1, 1)}{entropy};
    #Symmetry correction
    $entropy -= 1.4;
    my $r = 1.987; #molar gas constant
    my $tm = ($enthalpy * 1000 / ($entropy + ($r * log($oligo_conc))) - 
273.15 + (12* (log($salt_conc)/log(10))));
    return $tm;
}

sub format_printable_sequence
{
	my $seq = shift;


	$seq = join("\n", split /(.{80})/, $seq);
	$seq =~ s/\n\n/\n/g;
	return $seq;

}

sub get_feat
{
	my ($start, $end, $strand, $primary, $source_tag, $display_name, $label, $vntifkey) = @_;

	my $seqfeat = Bio::SeqFeature::Generic->new(
						-start=>	$start,
						-end=>		$end,
						-strand=>	$strand,
						-primary=>	$primary,
						-source_tag=>	$source_tag,
						-display_name=>	$display_name,
						-tag=>		{ 	label=>$label,
									vntifkey=>$vntifkey
								}
		);
	return $seqfeat;

}

sub create_gb
{
	my $name = shift;
	my $annot = shift;
	my $type = shift;
	my $seq_add = shift;
	my $fiveclipsize = shift;

	# Change this TODO
	my $infile = "/biodb/gateway/$type.gb";

	if(! -e $infile)
	{
		return '';
	}
	my $in = Bio::SeqIO->new(-file=>$infile, -format=>'genbank') or return;
	my $outseq_string = '';
	my $iostring = IO::String->new($outseq_string);
	my $out = Bio::SeqIO->new(-fh=>$iostring, -format=>'genbank');

	my $seq = $in->next_seq();

	my $seq_end = $seq->length();

	$seq->seq($seq->seq() . $seq_add);
	my $new_seqend = $seq->length();
	my @feats = ();

	push(@feats, get_feat($seq_end+1, $new_seqend, 0, 'misc_feature', 'vnti', 'TOPO insert', 'TOPO insert', 21) );
	push(@feats, get_feat($seq_end+1,  $seq_end + $fiveclipsize, 0, "5'clip", 'vnti', "gene 5' Clip", "gene 5' Clip", 51) );
	push(@feats, get_feat($seq_end + $fiveclipsize + 1, $new_seqend, 0, "CDS", 'vnti', "ORFid ", "Orfid $name - $annot", 4) );

	$seq->add_SeqFeature(@feats);
	$seq->display_id("pENTR/$type/ORF_" . $name . "-" . $fiveclipsize . "u");

	$out->write_seq($seq);
	return $outseq_string;

}




sub revcomp
{
	my $sequence = shift;
	$sequence =~ tr/ATGCatgcnN/TACGtacgnN/;
	$sequence = reverse($sequence);
	return $sequence;

}


