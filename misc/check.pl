#!/usr/bin/perl
  
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
 
use strict;
 
my $sequence_file = "query.fasta";
my $gene_file = "mcue.fasta";
 
my $gff_outfile = "out.gff";
 
 
# Create overlap array of hases
my $overlap;
my $overlap_bloodstream;
my $overlap_procyclic;
my $sequence_fasta_hash; 
 
my $genes =  Bio::SeqIO->new(   '-file'         => $gene_file,
                                '-format'       => "fasta");
 
my $sequences = Bio::SeqIO->new('-file'         => $sequence_file,
                                '-format'       => "fasta");
 
 
 
# Create blast databases for each sequence under the db directory
 
my $gene_hash;
open(GFF, ">", $gff_outfile);
 
open(GENEALL, ">", "db/all");
 
while(my $gene = $genes->next_seq)
{
        open(GENE, ">", "db/" . $gene->display_id);
 
        print GENE ">" . $gene->display_id . "\n" . $gene->seq();
        print GENEALL ">" . $gene->display_id . "\n" . $gene->seq() . "\n";
        close(GENE);
        system("cd db;formatdb -t " . $gene->display_id . " -i " . $gene->display_id . " -p F");
        $gene_hash->{$gene->display_id} = $gene;
        print GFF join("\t", $gene->display_id, "pcr", "gene", "1", $gene->length, ".", "+", "0", "Sequence " . $gene->display_id) . "\n";
        # initialize overlap array
        for(0..$gene->length)
        {
                $overlap->{$gene->display_id}[$_] = 0;
        }
 
 
}
 
system("cd db;formatdb -t all -i all -p F");
 
my $un_array;
my $ed_array;
 
my @all_params = (  'program' => 'blastn',
                    'F'=>'F',
                   'database' => "db/all",
                   "E" =>1e-3
);
 
my $all_factory = Bio::Tools::Run::StandAloneBlast->new(@all_params);
 
while(my $seq = $sequences->next_seq)
{
	$sequence_fasta_hash->{$seq->display_id} = $seq;

        if($seq->length < 10)
        {
                next;
        }
        print GFF join("\t", $seq->display_id, "pcr", "est", "1", $seq->length, ".", "+", "0", "Sequence " . $seq->display_id) . "\n";
        my $gff_type = 'dt_primed';
        if($seq->display_id() =~ /BS/)
        {
                $gff_type = 'specific_bloodstream';
        } elsif($seq->display_id() =~ /PC/)
        {
                $gff_type = 'specific_procyclic';
        }
        # First blast if vs the all db to find out which gene it is closest to
        my $all_report = $all_factory->blastall($seq);
        my $all_result = $all_report->next_result;
 
        my $hit = $all_result->next_hit;
        my $gene_name;
 
        if($hit)
        {
                my $hsp = $hit->next_hsp;
                my $unedited;
                my $edited;
                my $single_gene = 0;
                if($hsp->evalue < 1e-3)
                {
                } else
                {
                        next;
                }
 
                my $gene_hit =  $hit->name;
 
                if($gene_hit =~ /ed$/)
                {
                        ($gene_name) = $gene_hit =~ /(.+)ed$/;
                        $unedited = $gene_name . "un";
                        $edited = $gene_name . "ed";
                } elsif($gene_hit =~ /un$/)
                {
                        ($gene_name) = $gene_hit =~ /(.+)un$/;
                        $unedited = $gene_name . "un";
                        $edited = $gene_name . "ed";
                }else
                {
                        $gene_name = $gene_hit;
                        $unedited = $gene_name;
                        $edited = $gene_name;
                        $single_gene = 1;
                }
 
 
                my $unedit_seq = $gene_hash->{$unedited};
                my $edit_seq = $gene_hash->{$edited};
 
 
                my @un_params = (  'program' => 'blastn',
                                        'F'=>'F',
                                'database' => "db/" . $unedited
                );
 
                my @ed_params = (  'program' => 'blastn',
                                        'F'=>'F',
                                'database' => "db/" . $edited
                );
 
                my $un_factory = Bio::Tools::Run::StandAloneBlast->new(@un_params);
                my $ed_factory = Bio::Tools::Run::StandAloneBlast->new(@ed_params);
                                                                                                                                                                                                                                                    
                my $max_unedit_size = 0;
                my $max_edit_size = 0;
                my $max_edit_start = 0;
                my $max_edit_stop = 0;
                my $max_unedit_start = 0;
                my $max_unedit_stop = 0;
                my $un_hash;
                my $ed_hash;
 
 
                if($gene_hash->{$unedited})
                {
                        # Find area's
                        my $un_report = $un_factory->blastall($seq);
                        my $un_result = $un_report->next_result;
 
                        while(my $hit = $un_result->next_hit)
                        {
                                while(my $hsp = $hit->next_hsp)
                                {
 
                                        if( ($hsp->hit->end - $hsp->hit->start) > $max_unedit_size)
                                        {
                                                $max_unedit_size = $hsp->hit->end - $hsp->hit->start;
                                                $max_unedit_start = $hsp->hit->start;
                                                $max_unedit_stop = $hsp->hit->end;
                                        }
 
                                }
                        }
                        $un_hash->{size} = $max_unedit_size;
                        $un_hash->{start} = $max_unedit_start;
                        $un_hash->{stop} = $max_unedit_stop;
                        $un_hash->{name} = $seq->display_id;
 
                        print GFF join("\t", $unedit_seq->display_id, "HSP", $gff_type, $un_hash->{start}, $un_hash->{stop}, ".", ".", ".", "Sequence " . $un_hash->{name}) . "\n";
                        for($un_hash->{start}..$un_hash->{stop})
                        {
                                $overlap->{$unedit_seq->display_id}[$_]++;
                        }
                }
 
                if($gene_hash->{$edited} && !($single_gene))
                {
                        my $ed_report = $ed_factory->blastall($seq);
                        my $ed_result = $ed_report->next_result;
 
                        while(my $hit = $ed_result->next_hit)
                        {
                                while(my $hsp = $hit->next_hsp)
                                {
                                                                                                                                                                                                                                                    
               
                                        if( ($hsp->hit->end - $hsp->hit->start) > $max_edit_size)
                                        {
                                                $max_edit_size = $hsp->hit->end - $hsp->hit->start;
                                                $max_edit_start = $hsp->hit->start;
                                                $max_edit_stop = $hsp->hit->end;
                                        }
                                }
                        }
                        $ed_hash->{size} = $max_edit_size;
                        $ed_hash->{start} = $max_edit_start;
                        $ed_hash->{stop} = $max_edit_stop;
                        $ed_hash->{name} = $seq->display_id;
 
                        print GFF join("\t", $edit_seq->display_id, "HSP", $gff_type, $ed_hash->{start}, $ed_hash->{stop}, ".", ".", ".", "Sequence " . $ed_hash->{name}) . "\n";
                        for($ed_hash->{start}..$ed_hash->{stop})
                        {
                                $overlap->{$edit_seq->display_id}->[$_]++;
                        }
                }
 
        } # End If hit
} # End while seq
 
 
# Now print out the overlap for the GFF
 
while(my ($gene_name, $overlap_array) = each(%$overlap))
{
        my $count = 0;
        foreach my $overlap_score (@$overlap_array)
        {
                if($count == 0)
                {
                } else
                {
                        print GFF join("\t", $gene_name, "overlap", "overlap", $count, $count, $overlap_score, ".", ".") . "\n";
                }
                $count++;
        }
}

