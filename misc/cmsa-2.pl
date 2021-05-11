#!/usr/bin/perl -w


use strict;
use Mbl;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
use IO::String;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::AlignIO;
use Bio::DB::GenPept;
use Bio::Tools::BPpsilite;

# insert into phylo_pipe (id, idname, sequence_type_id, translate, sequence, type) SELECT NULL, orfid, 2, 'Y', sequence, 'psiblast' from orfs where delete_fg = 'N';

#$ENV{'BLASTDB'} = '/blastdb';

my $min_eval = 1e-7;


my $mbl = Mbl::new(undef, 'giardia');
my $dbh = $mbl->dbh();
my $fastadir = "/bioware/fasta34/";
my $blastdb_dir = "/blastdb/";
my $orfdb = "giardia12_orfs_aa";
my $not_include_taxa = "184922";

my $sth = $dbh->prepare("select id as id, idname , sequence, translate from phylo_pipe where type='psiblast' LIMIT 1");
my $del_pp = $dbh->prepare("delete from phylo_pipe where id = ?");
my $ins_orthopara = $dbh->prepare("insert into orthopara (idname, db_id, db_xref) VALUES (?, ?, ?)");
my $del_orthopara = $dbh->prepare("delete from orthopara where idname = ?");
my $ins_ma = $dbh->prepare("insert into ma (id, idname, ma, type, description) VALUES (?, ?, ?, ?, ?)");
my $del_ma = $dbh->prepare("delete from ma where idname = ? and type = ?");
my $ins_tree = $dbh->prepare("insert into tree (ma_id, idname, tree, type, description) VALUES (?, ?, ?, ?, ?)");
my $del_tree = $dbh->prepare("delete from tree where idname = ? AND type = ?");


my $keep_going = 1;
while($keep_going)
{
	$sth->execute();
	if($sth->rows > 0)
	{
		my $row = $sth->fetchrow_hashref;
		my $visited_orfs = 0;

		my $orfid = $row->{idname};
		$del_pp->execute($row->{id});

		mkdir($orfid);

		my @orf_params = ( 'program'=> 'blastp',
				'database'=> "$orfdb",
		                'e'=>$min_eval);

		my @nr_params = ( 'program'=> 'blastp',
		                'database'=> $blastdb_dir . "nr",
				'I'=> 'T',
				'e'=>$min_eval);


		my @orf_psi_params = (
		                'database'=> $orfdb,
		                'j'=> '3',
				'e'=>'1e-3',
		                'h'=>$min_eval);

		my @psi_params = ( 
		                'database'=> $blastdb_dir . "nr",
		                'j'=> '2',
		                'I'=> 'T',
				'R'=>"$orfid/$orfid.cp.1",
		                'h'=>$min_eval);

		my $orf_factory = Bio::Tools::Run::StandAloneBlast->new(@orf_psi_params);
	
		my $nr_factory =  Bio::Tools::Run::StandAloneBlast->new(@nr_params);

		my $protein_input = orfid_to_seq_aa($orfid);
		my @tmp;
		my @seq_array_aa;

		push(@tmp, $protein_input);
		seq_to_fasta(\@tmp, $orfid . '/' . $orfid . ".fas");
	
		push(@seq_array_aa, $protein_input);

		# Blast against own orfs

		my $blast_report = $orf_factory->blastpgp($protein_input);

		my $this_result = $blast_report->next_result;

		while(my $hit = $this_result->next_hit())
		{
			my $seq = orfid_to_seq_aa($hit->name);
			if($seq->display_id ne  "gn1|gdb|" . $orfid )
			{
	                        my $result_eval = compare_prss($protein_input, $seq, 1000, 8080);
	                        if($result_eval < ($min_eval*1000))
	                        {
	                                push(@seq_array_aa, $seq);
	                        }
			}
		}

		# Now create am temporary blast database of these orfs
		seq_to_fasta(\@seq_array_aa, $orfid . '/' . $orfid);
	
		system("formatdb -t $orfid -i $orfid/$orfid -p T ");

		system("blastpgp -d $orfid/$orfid -e 10000 -j 2 -C $orfid/$orfid.cp.1 -Q $orfid/$orfid.cp.1.txt -i $orfid/$orfid.fas -o $orfid/$orfid.initial.psi");


		my $psi_factory =  Bio::Tools::Run::StandAloneBlast->new(@psi_params);
		my $psi_report = $psi_factory->blastpgp($protein_input);

		my $psi_result = $psi_report->next_result;
		my $num_dbentries = $psi_result->database_entries();

		while(my $hit = $psi_result->next_hit())
		{
			my ($gi_num) = $hit->name =~ /gi\|(\d+)/;
			my $taxid = $mbl->get_taxid_from_gi($gi_num);
			if($taxid)
			{
				if($taxid eq $not_include_taxa)
				{
					# Do nothing
				} else
				{
			        	my $seq = genbank_to_seq_aa($gi_num);
				        push(@seq_array_aa, $seq);
				}
			}
		}

		my @new_seq_array_aa;

		my $prot_seq = shift(@seq_array_aa);
		push(@new_seq_array_aa, $prot_seq);

		foreach my $seq (@seq_array_aa)
		{
			my $result_eval = compare_prss($prot_seq, $seq, 1000, $num_dbentries);
			if($result_eval < ($min_eval*1000))
			{
				push(@new_seq_array_aa, $seq);
			}
		}

		@seq_array_aa = @new_seq_array_aa;

		my $array_ref = \@seq_array_aa;
		# remove duplicates
		$array_ref = remove_duplicate_seq($array_ref);

		$del_orthopara->execute($orfid);
		if(scalar @seq_array_aa > 1)
		{
			foreach my $seq (@seq_array_aa)
			{
				my $id;
				my $db;
				if($seq->display_id() =~ /^gi/)
				{
					$db = 2;
					($id) = $seq->display_id() =~ /^gi\|(\d+)/;
				} elsif($seq->display_id() =~ /^gn1/)
				{
					$db = 7;
					($id) = $seq->display_id() =~ /^gn1\|gdb\|(\d+)/;
				} else
				{
					#
				}
				$ins_orthopara->execute($orfid, $db, $id);
			}
			# Create alignment
	
			my @align_params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'OUTPUT'=>'GCG');
			my $align_factory = Bio::Tools::Run::Alignment::Clustalw->new(@align_params);
			$align_factory->quiet(1);

			my $aln = $align_factory->align($array_ref);

			my $align_str = create_alignment_string($aln, 'nexus');
			create_alignment_file($aln, 'clustalw', $orfid . "/" . $orfid . ".aln");
			$del_ma->execute($orfid, 'nexus');
			$ins_ma->execute($row->{id}, $orfid, $align_str, 'nexus', 'Automaticaly Generated using psiblast vs own genome, then psiblast vs nr and alignment created with ClustalW');

			# Create tree now
			system("cd $orfid;clustalw -tree -outputtree=nj -infile=" . $orfid . ".aln");
			my $njtree = get_string_from_file($orfid . '/' . $orfid . ".ph");
			$del_tree->execute($orfid, 'nj');
			$ins_tree->execute($row->{id}, $orfid, $njtree, 'nj', 'Automaticaly Generated using psiblast vs own genome, then psiblast vs nr and alignment and tree created with ClustalW');

		}

		system("rm -rf $orfid");

	} else
	{
		$keep_going = 0;

	}# END ALL

} # END WHILE KEEP GOING





sub orfid_to_seq_nt
{
	my $orfid = shift;
	my $q = $dbh->prepare('select orfid as orfid, sequence, start, stop, direction from orfs where orfid = ?');
	$q->execute($orfid);
	my $row = $q->fetchrow_hashref;
	my $retseq = Bio::Seq->new( -display_id=>"gn1|gdb|" .$orfid , -seq=>$row->{sequence} );
	return $retseq;

}

sub orfid_to_seq_aa
{
        my $orfid = shift;
        my $q = $dbh->prepare('select orfid as orfid, sequence, start, stop, direction from orfs where orfid = ?');
        $q->execute($orfid);
        my $row = $q->fetchrow_hashref;
	my $retseq = Bio::Seq->new( -display_id=>"gn1|gdb|" . $orfid , -seq=>$row->{sequence} );

	$retseq = $retseq->translate();
        return $retseq;
                                                                                                                        
}

sub genbank_to_seq_aa
{
	my $gi = shift;
	`fastacmd -d /blastdb/nr -s $gi > .$gi`;
	my $seqio = Bio::SeqIO->new(-file =>".$gi",-format=>"fasta");
	my $seq = $seqio->next_seq();
	my $retseq =  Bio::Seq->new( -display_id=>'gi|'. $gi , -seq=>$seq->seq);
	`rm -f .$gi`;
	return $retseq;
}

sub seq_to_fasta
{
	my $seqs = shift;
	my $filename = shift;
	my @seq_array = @{$seqs};

	open(FASTA, ">", $filename);

	foreach my $seq (@seq_array)
	{
		if(!$seq->seq())
		{
			print "NO SEQUENCE FOUND FOR " . $seq->display_id . "\n";
		} else
		{	
			print FASTA ">" . $seq->display_id . "\n";
			print FASTA $seq->seq() . "\n";
		}
	}
	close(FASTA);
}

sub compare_prss
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $iterations = shift;;
	my $num_entries = shift;

	my $seq1_name = int(rand(100000));
	my $seq2_name = int(rand(100000));

	open(FASTA1, ">", "/tmp/.$seq1_name.prss");
        open(FASTA2, ">", "/tmp/.$seq2_name.prss");

	print FASTA1 ">$seq1_name\n";
	print FASTA1 $seq1->seq();
        print FASTA2 ">$seq2_name\n";
        print FASTA2 $seq2->seq();
	close(FASTA1);
	close(FASTA2);

	my $result = `/opt/BioBrew/Fasta/3.4/bin/prss34 -s BL62 -f 11 /tmp/.$seq1_name.prss /tmp/.$seq2_name.prss $iterations`;

	my ($eval) = $result =~ /unshuffled.+\<(.+)/;
	$eval =~ s/ //gi;

	system('rm -f /tmp/.' . $seq1_name . '.prss');
        system('rm -f /tmp/.' . $seq2_name . '.prss');

	return ($eval * ($num_entries+$iterations));
	

}

sub remove_duplicate_seq
{
	my $seq_array = shift;
	my %seq_hash;
	foreach my $seq (@{$seq_array})
	{
		$seq_hash{$seq->display_id} = $seq;
	}
	my $ret_array;
	while(my ($key, $val) = each %seq_hash)
	{
		push(@{$ret_array}, $val);
	}
	return $ret_array;

}

sub create_alignment_file
{
	my $aln = shift;
	my $format = shift;
	my $filename = shift;

	my $str;
        my $out = IO::String->new(\$str);
        my $ioout = Bio::AlignIO->new(-format=> $format, -fh => $out );
	$ioout->write_aln($aln);
        open (ALIGN, ">" . $filename);
	print ALIGN $str;
        close(ALIGN);
	return 1;

}

sub create_alignment_string
{
        my $aln = shift;
        my $format = shift;

        my $str;
        my $out = IO::String->new(\$str);
        my $ioout = Bio::AlignIO->new(-format=> $format, -fh => $out );
        $ioout->write_aln($aln);
	return $str;		
}

sub find_local_orfs
{
	my $seq_array = shift;
	my $new_array;
	foreach my $seq (@{$seq_array})
	{
		if($seq->display_id =~ /gdb/)
		{
			my ($orfid) = $seq->display_id	=~ /gn1\|gdb\|(\d+)/;
			push(@{$new_array}, $orfid);
		}
	}
	return $new_array;
}

sub get_string_from_file
{
	my $filename = shift;
	open(FILE, $filename);
	my $str = '';

	while(<FILE>)
	{
		$str .= $_;
	}
	return $str;
}
