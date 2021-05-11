#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Restriction::Analysis;
use Bio::Restriction::IO;
use Data::Dumper;


# Had to fix /usr/local/share/perl5/Bio/Restriction/IO/withrefm.pm to call -vendors instead of -vendor which does not work
use strict;

my $rebase = Bio::Restriction::IO->new(-file=>'./withrefm.310', -format=>'withrefm', -verbose=>1);
my $rebase_collection = $rebase->read();
my $commercial_collection = Bio::Restriction::EnzymeCollection->new( -empty => 1 );
$commercial_collection->enzymes(grep {has_vendor($_)} $rebase_collection->each_enzyme() );

# This will print all commercial enzymes
#map { print $_->name . "\n";} $commercial_collection->each_enzyme();
#exit;
#
foreach my $fasta_file(@ARGV)
{
	my $seqs = Bio::SeqIO->new(-file=>$fasta_file);

	while(my $seq = $seqs->next_seq)
	{	

		my $ra = Bio::Restriction::Analysis->new(-enzymes=>$rebase_collection, -seq=>$seq);
		my $cutters = $ra->cutters();
	
		map {
			
			my $enz = $rebase_collection->get_enzyme($_->name());
			
			if(has_vendor($enz))
			{
				# print $enz->name() . " " . join(" ", @{get_vendors($enz)});# Do nothing
				print $enz->name() . "\n";
				#print $seq->display_id . "\t" . $enz->name() . "\t" . join(",", $ra->sizes($enz)) . "\n";
				# print $seq->display_id . "\t" . $enz->name() . "\t" . join(",", $ra->positions($enz->name())) . "\n";
			} else
			{
				# print " " . $enz->name(); # Do nothing
			}

		    } $cutters->each_enzyme;
	}
}



sub has_vendor
{
	my $enz = shift;
	my $specific_vendor = shift;

	if(scalar @{get_vendors($enz)} > 0)
	{
		if(defined($specific_vendor))
		{
			if($specific_vendor ~~ @{get_vendors($enz)})
			{
				return 1;
			} else
			{
				return 0;
			}
		}
		return 1;
	} else
	{
		return 0;
	}
}

sub get_vendors
{
	my $enz = shift;
	my @vendors = ();
	if(ref($enz->{_vendors}) eq 'ARRAY')
	{
		push(@vendors,  @{$enz->{_vendors}->[0]});
		return \@vendors;
	} else
	{
		return \();
	}
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
           (A parameter should be added to allow the oligo concentration 
to be set.)
 Notes   : Calculation of Tm as per Allawi et. al Biochemistry 1997 
36:10581-10594.  Also see
           documentation at http://biotools.idtdna.com/analyzer/ as they 
use this formula and
           have a couple nice help pages.  These Tm values will be about 
are about 0.5-3 degrees
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


