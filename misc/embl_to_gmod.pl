#!/usr/bin/perl

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Getopt::Long;
use File::Temp qw/ tempdir tempfile/;

use Mbl;
use DBI;


use strict;


my $contig_file;
my $output_prefix;
my $verbose = 1;
my $gmoddb;
my $gff_file;

my $options = GetOptions(
				"contig_file=s", \$contig_file,
				"output_prefix=s", \$output_prefix,
				"organism_db=s", \$gmoddb,
				"gff_file=s", \$gff_file,
				"verbose!", \$verbose
			);




if(!defined($contig_file) || !defined($output_prefix) || !defined($gff_file) )
{
	print "
--contig_file
--organism_db
--output_prefix
--gff_file
--verbose

Example: program.pl --contig_file=t_cruzi_contigs.embl --output_prefix=t_cruzi

";
exit;
}


my $mbl = Mbl::new(undef, $gmoddb);
my $dbh = $mbl->dbh;
my $insert_database = 1;

# Set all database queries
my $insert_contig               = $dbh->prepare('insert into links (super_id, bases_in_super, contigs_in_super, ordinal_number, contig_length, gap_before_contig, gap_after_contig, contig_number, contig_start_super_base, modified_contig_start_base, modified_bases_in_super) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
my $insert_contig_seq           = $dbh->prepare('insert into contigs (contig_number, bases) values (?, ?)');
my $insert_read                 = $dbh->prepare('insert into reads (read_name) values (?)');
my $insert_into_read_assembly   = $dbh->prepare('insert into reads_assembly (read_name, read_len_untrim, first_base_of_trim, read_len_trim, contig_number, contig_length, trim_read_in_contig_start, trim_read_in_contig_stop, orientation) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)');
my $insert_reads_bases          = $dbh->prepare('insert into reads_bases (read_name, bases) VALUES (?, ?)');
my $insert_orf                  = $dbh->prepare('insert into orfs (orfid, orf_name, annotation, annotation_type, source, contig, start, stop, direction, delete_fg, delete_reason, sequence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
my $insert_annotation           = $dbh->prepare('insert into annotation (userid, orfid, update_dt, annotation, notes, delete_fg, blessed_fg, qualifier, with_from, aspect, object_type, evidence_code, private_fg) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
        

open(GFFFILE, ">", $gff_file);                                                                                                                                  

my $tempdir = tempdir(CLEANUP=>1);


my $next_scaffold_num = 1;
my $next_contig_num = 1;
my $scaffold_lookup_hash;
my $contig_lookup_hash;
my $scaffold_hash;
my $contig_hash;

# A scaffold will have the following (scaffold_id, scaffold_accession, scaffold_accession_version, contig_array, scaffold_size, scaffold_gapped_size)
# A scaffold_array will be an array of contigs with the following values (contig_id, contig_accession, contig_accession_version, complement, start, stop, gap_before, gap_after, ordinal_number, length, contig_start_scaffold, contig_start_scaffold_gapped)

my $scaffold_record;
my $scaffold_accession;
my $scaffold_accession_version;

open(CONTIGFILE, $contig_file);
my $contig_line;
my $contig_accession;
my $contig_accession_version;
my $contig_single_file = "";
while(<CONTIGFILE>)
{
        my $contig_line = $_;
        $contig_single_file .= $contig_line;
        if($contig_line =~ /^\/\//)
        {
                # print out the record to a new file
                open(CONTIG_TMP_FILE, ">" . $tempdir . "/" . $contig_accession . ".embl");
                print CONTIG_TMP_FILE $contig_single_file;
                close(CONTIG_TMP_FILE);
                # reset the variables
                $contig_accession = undef;
                $contig_accession_version = undef;
                $contig_single_file = "";
        } elsif($contig_line =~ /^ID/)
        {
                ($contig_accession) = $contig_line =~ /^ID\s+(\S+)\s*/;
        }
                                                                                                                                           
}
close(CONTIGFILE);





my $contigcheck_seqio = Bio::SeqIO->new(-file=>$contig_file, -format=>'embl');

while(my $seq = $contigcheck_seqio->next_seq)
{
	print $seq->id . "\n";
	if(defined($contig_hash->{$seq->display_id}))
	{
		# Do nothing
	} else
	{
		# Create a new scaffold with only this one contig
		my $scaffold_info;
		my $scaffold_num = $next_scaffold_num;
	        $scaffold_info->{scaffold_id} = $scaffold_num;
		$scaffold_info->{scaffold_accession} = $seq->display_id;
		$scaffold_info->{scaffold_accession_version} = 0;
	        $next_scaffold_num++;
		my $contig_array;
                my $contig_record;
		my $contig_num = $next_contig_num;
		$next_contig_num++;
                $contig_record->{contig_id} = $contig_num;
                $contig_record->{contig_accession} = $seq->display_id;
                $contig_record->{contig_accession_version} = 0;
                $contig_record->{complement} = 0;
                $contig_record->{start} = 1;
                $contig_record->{stop} = $seq->length;
                $contig_record->{length} = $seq->length;
                $contig_record->{gap_before} = 0;
                $contig_record->{gap_after} = 0;
                $contig_record->{ordinal_number} = 1;
                $contig_record->{contig_start_scaffold} = 0;
                $contig_record->{contig_start_scaffold_gapped} = 0;
        	$scaffold_info->{scaffold_size} = $seq->length;
	        $scaffold_info->{scaffold_gapped_size} = $seq->length;

                push(@$contig_array, $contig_record);
                $contig_hash->{$seq->display_id} = $scaffold_num;
		$scaffold_info->{contig_array} = $contig_array;
		$scaffold_hash->{$seq->display_id} = $scaffold_info;
	}
}


my @sorted_scaffolds = sort { $scaffold_hash->{$a}->{scaffold_id} <=> $scaffold_hash->{$b}->{scaffold_id} } keys %$scaffold_hash;

foreach my $scaffold_acc (@sorted_scaffolds)
{
	my $scaffold = $scaffold_hash->{$scaffold_acc};
	if($verbose)
	{
		print join("\t", "IMPORTING:", $scaffold_acc, $scaffold->{scaffold_id}) . "\n";
	}
	
	insert_scaffold($scaffold);

}


if($verbose)
{
	print "Complete!\n";
}



sub contig_lookup
{
	my $accession_number = shift;

	my $contig_file_name = $tempdir . "/" . $accession_number . ".embl";

	my $contig_seqio = Bio::SeqIO->new(-file=>$contig_file_name, -format=>'embl');
	my $contig_seq = $contig_seqio->next_seq;

	return $contig_seq;
	
}

sub insert_scaffold
{
	my $scaffold = shift;


	my $num_contigs = scalar @{$scaffold->{contig_array}};
	foreach my $contig( @{$scaffold->{contig_array}})
	{
		# First insert into the links table
		if($insert_database)
		{
			$insert_contig->execute($scaffold->{scaffold_id}, $scaffold->{scaffold_size}, $num_contigs, $contig->{ordinal_number}, $contig->{length}, $contig->{gap_before}, $contig->{gap_after}, $contig->{contig_id}, $contig->{contig_start_scaffold}, $contig->{contig_start_scaffold_gapped}, $scaffold->{scaffold_gapped_size});
		}
		# Now insert into the contigs table (sequence)

		my $seq = contig_lookup($contig->{contig_accession});

		if($insert_database)
		{
			# If it is complemented, then complement it here
			if($contig->{complement})
			{
				$insert_contig_seq->execute("contig_" . $contig->{contig_id}, $seq->revcom()->seq());
			} else
			{
				$insert_contig_seq->execute("contig_" . $contig->{contig_id}, $seq->seq());
			}
		}
		# Now insert it as a read in the original orientation
		my $direction = "+";
		if($contig->{complement} > 0)
		{
			$direction = "-";
		}
		if($insert_database)
		{
			$insert_read->execute($contig->{contig_accession});
			$insert_into_read_assembly->execute($contig->{contig_accession}, $contig->{length}, 1, $contig->{length}, $contig->{contig_id}, $contig->{length}, 1, $contig->{stop}, $direction);
			$insert_reads_bases->execute($contig->{contig_accession}, $seq->seq());
		}

		# Now we have the read and the contigs in, so we just need to import the orf calls

		insert_orf_information($seq, $contig->{contig_id}, $contig->{complement}, $contig->{length});

		# TODO: output a embl2gff file of all of the annotations on this contig

		emblgff_contig($seq, $contig->{length}, $contig->{contig_id}, $scaffold->{scaffold_id}, $contig->{complement});
	}

}


sub emblgff_contig
{
	my $seq = shift;
	my $contig_length = shift;
	my $contig_id = shift;
	my $scaffold_id = shift;
	my $complement = shift;


	foreach my $feat ($seq->get_SeqFeatures)
	{
		my $gff_string = $feat->gff_string();
		my $contig_string = 'contig_' . $contig_id;
		$gff_string =~ s/^SEQ/$contig_string/;
		my @gff = split("\t", $gff_string);
		# Check if this is the complement, and if it is, reverse the direction
		if($complement)
		{
			my $feat_start = $gff[3];
			my $feat_stop = $gff[4];
			my $dir = $gff[6];
			($feat_start, $feat_stop, $dir) = get_coordinates($feat_start, $feat_stop, $dir, $complement, $contig_length);
			$gff[3] = $feat_start;
			$gff[4] = $feat_stop;
			$gff[6] = $dir;
		}
		# Now we have to 
		print GFFFILE join("\t", @gff) . "\n";
		
	}
}

sub insert_orf_information
{
	
	my $seq = shift;
	my $contig_id = shift;
	my $complement = shift;
	my $contig_length = shift;

	foreach my $feat ($seq->top_SeqFeatures)
	{
		my $dir = "+";
		if($feat->strand == -1)
		{
			$dir = '-';
		} elsif($feat->strand == 0)
		{
			$dir = '.';
		}
		my $feat_start = $feat->start;
		my $feat_stop = $feat->end;
		my $feat_sequence = $feat->seq()->seq();
		if($feat->start > $feat->end)
		{
			my $t = $feat_start;
			$feat_start = $feat_stop;
			$feat_stop = $t;
		}
                # Is this a gene
                my $feat_hash;
                my $note_text = '';
                my $gene = '';
                my $product = '';
		my $codon_start = 1;
                foreach my $tag ($feat->get_all_tags() )
                {
                                                                                                                                          
                        if($tag eq 'gene')
                        {
                                $gene = join(". ", $feat->get_tag_values($tag));
                                                                                                                                          
                        } elsif($tag eq 'product')
                        {
                                $product = join(". ", $feat->get_tag_values($tag));
                        }elsif($tag eq 'codon_start')
			{
				$codon_start = $feat->get_tag_values($tag);
				if($codon_start > 1)
				{
					if($dir eq "+")
					{
						$feat_start = $feat_start + $codon_start - 1;
					} elsif($dir eq "-")
					{
						$feat_stop = $feat_stop - $codon_start + 1;
					}
					my $feat_seq_obj = $feat->seq();
					$feat_sequence = $feat_seq_obj->subseq($codon_start, $feat_seq_obj->length());
				}
			} else
                        {
                                $note_text .= $tag . ':' . join(". \n", $feat->get_tag_values($tag)) . "\n";
                        }
                                                                                                                                          
                }

                # If this was a gene, then lets create an orf for it.
                if( $feat->primary_tag eq 'CDS' )
                {
                        if($insert_database)
                        {
                                my $orfid = get_max_orfid()+1;
                                if($verbose)
                                {
                                        print "##CREATING ORF $orfid##\n";
                                }
                                my $gene_name = $gene;
                                if($gene_name eq '')
                                {
                                        $gene_name = $product;
                                } elsif($product ne '')
                                {
                                        $gene_name = $gene . ' - ' . $product;
                                } else
                                {
                                        $gene_name = 'unnamed';
                                }
				if($complement)
				{
					if($verbose)
					{
						print "COMPLEMENT\n";
						print join("\t", $feat_start, $feat_stop, $dir, $complement, $contig_length);
					}
					($feat_start, $feat_stop, $dir) = get_coordinates($feat_start, $feat_stop, $dir, $complement, $contig_length);
					if($verbose)
					{
						print join("\t", $feat_start, $feat_stop, $dir, $complement, $contig_length);
					}
				
				}
				if($insert_database)
				{
	                                $insert_orf->execute($orfid, undef, undef, undef, 'GENBANK', 'contig_' . $contig_id, $feat_start, $feat_stop, $dir, 'N', undef, $feat_sequence);
	                                $insert_annotation->execute(1, $orfid, undef, $gene_name, $note_text, 'N', 'Y', undef, undef, undef, undef, undef, 'N');
				}
                        }
                                                                                                                                          
                }


	}
}
sub get_coordinates
{
	my $start = shift;
	my $stop = shift;
	my $direction = shift;
	my $complement = shift;
	my $length = shift;

	if(!$complement)
	{
		return ($start, $stop, $direction);
	} else
	{
		my $new_stop = $length - $start +1;
		my $new_start = $length - $stop +1;
		my $new_direction = "+";
		if($direction eq "+")
		{
			$new_direction = "-";
		}
		return( $new_start, $new_stop, $new_direction);
	}
}

sub get_max_orfid
{
        my $sth = $dbh->prepare("select max(orfid) as max_id from orfs");
        $sth->execute();
        return $sth->fetchrow_hashref->{max_id};
}

