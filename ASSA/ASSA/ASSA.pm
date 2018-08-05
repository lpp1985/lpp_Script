package ASSA::ASSA;

use strict;
use warnings;

# $Id: ASSA.pm 2581 2016-09-13 14:20:59Z antonov $

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
# 
# !!! ALL COORDINATES ARE 1-BASED !!!
#
# $self
#   |
#   |
#   |--->{unique_id}     --  next unique ID
#   |--->{bdb_path}      --  path to the BLASTn db
#   |--->{tmp_dir}
#   |--->{verbose}
#   |
#   |
#   |--->{q_trx_hash}
#   |--->{h_trx_hash}
#   |          |
#   |          |--->{$trx_id}
#   |                   |
#   |                   |--->{TRX_ID}
#   |                   |--->{GENE_ID}
#   |                   |--->{SEQ}
#   |                   |
#   |                   |--->{SEQ_ALU}
#   |                   |--->{NUM_ALU}
#   |                   |--->{CIS_GENES}  --  [hash with GENE_IDs] added after remove_cis_gene_pairs() was called
#   |
#   |--->{q_gene_hash}
#   |--->{h_gene_hash}
#   |           |
#   |           |--->{$gene_id}   --   array of trx hashes
#   |
#   |
#   |
#   |---{blast_sites}        --  will be generated after the run_blast()
#   |         |
#   |         |--->{$q_gene_id}
#   |                    |
#   |                    |--->{$h_gene_id}   --  hash of sites for the best TRX pair
#   |                               |
#   |                               |--->{$site_id}
#   |                                          |
#   |                                          |--->{ID}          --  unique site ID
#   |                                          |--->{Q_GENE_ID}
#   |                                          |--->{H_GENE_ID}
#   |                                          |--->{Q_TRX_ID}
#   |                                          |--->{H_TRX_ID}
#   |                                          |--->{Q_TRX_LEN}
#   |                                          |--->{H_TRX_LEN}
#   |                                          |--->{Q_START}
#   |                                          |--->{Q_END}
#   |                                          |--->{H_START}
#   |                                          |--->{H_END}
#   |                                          |--->{ALI_LEN}
#   |                                          |--->{GC_CONTENT}  --  [0..100]
#   |                                          |--->{IDENTITY}    --  [70..100]
#   |                                          |--->{GAPS_OPENINGS}
#   |                                          |--->{MISMATCHES}
#   |                                          |--->{SCORE}
#   |                                          |--->{EVALUE}
#   |
#   |
#   |
#   |---{bifold_sites}        --  will be generated after the run_bifold()
#   |            |
#   |            |--->{$q_gene_id}
#   |                       |
#   |                       |--->{$h_gene_id}
#   |                                  |
#   |                                  |--->{$site_id}
#   |                                             |
#   |                                             |--->{ID}              --  
#   |                                             |--->{BLAST_SITE_ID}   --  
#   |                                             |
#   |                                             |--->{Q_GENE_ID}  -- ?
#   |                                             |--->{H_GENE_ID}  -- ?
#   |                                             |--->{Q_TRX_ID}
#   |                                             |--->{H_TRX_ID}
#   |                                             |--->{Q_SEQ}
#   |                                             |--->{Q_START}
#   |                                             |--->{Q_END}
#   |                                             |--->{H_SEQ}
#   |                                             |--->{H_START}
#   |                                             |--->{H_END}
#   |                                             |
#   |                                             |--->{DPL_LEN}
#   |                                             |--->{DPL_LEN_PVALUE}
#   |                                             |--->{DDG}
#   |                                             |--->{DDG_PVALUE}
#   |                                             |--->{DPL_ARR}
#   |
#   |
#   |---{coexpression}        --  will be generated after the filter_gene_pairs_by_coexpression()
#   |            |
#   |            |--->{$q_gene_id}
#   |                       |
#   |                       |--->{$h_gene_id}
#   |                                  |
#   |                                  |--->{PCC}         --  pearson correlation coefficient
#   |                                  |--->{PCC_PVALUE}
#   |                                  |--->{PCC_NUM}     --  number of expression values that were used to compute the PCC
#   |
#

use Data::Dumper;
use Carp qw(confess); 
use IPC::Open2;
use POSIX qw(:sys_wait_h);
use Cwd qw(abs_path);

use ASSA::Lib qw(write_to_file ah2a fasta2hash intersect_digits read_table_wo_header unique sum min max log10);
use ASSA::Gauss;
use ASSA::Ct_formatManager;

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw(
	compute_ddG_pvalue
	compute_dpl_len_pvalue
);

###
# CONSTANTS


###
# CONSTRUCTOR
#
# %opts
#   |
#   |--->{q_trx_hash}
#   |--->{h_trx_hash}
#   |          |
#   |          |--->{$trx_id}
#   |                   |
#   |                   |--->{TRX_ID}
#   |                   |--->{GENE_ID}
#   |                   |--->{SEQ}
# 
#
sub new
{
	my $class = shift;
	my(%o) = @_;
	
	die "Either q_trx_hash or q_fn option is required: ".Dumper(\%o) if !$o{q_trx_hash} && !$o{q_fn};
	die "Either h_trx_hash or h_fn option is required: ".Dumper(\%o) if !$o{h_trx_hash} && !$o{h_fn};
	die "tmp_dir is required: ".Dumper(\%o)                          if !$o{tmp_dir};
	
	mkdir($o{tmp_dir}, 0777) if !-d $o{tmp_dir};
	mkdir($o{out_dir}      ) if     $o{out_dir};
	
	my $self = bless {
		q_trx_hash  => $o{q_trx_hash} || fasta2hash($o{q_fn}),
		h_trx_hash  => $o{h_trx_hash} || fasta2hash($o{h_fn}),
		q_gene_hash => {},
		h_gene_hash => {},
		gene_pairs  => {},
		unique_id   => 1,
		cis_seed    => 30,
		num_structs => 20,      # num bifold structs
		num_threads => $o{num_threads} || 1,
#		delim       => '_',
		bdb_path    => "$o{tmp_dir}/_target_blast_db",
		out_dir     => $o{out_dir},
		tmp_dir     => $o{tmp_dir},
		verbose     => $o{verbose},
	}, $class;
	
	# verification
	foreach my $trx ( values(%{$self->{q_trx_hash}}) , values(%{$self->{h_trx_hash}}) )
	{
		warn "The sequence should consist of ACGT chars only: ".Dumper($trx) if $trx->{SEQ} =~ /[^ACGT]/i;
	}
	
	add_alu_info_to_trx_hash($self->{q_trx_hash}, $o{q_alu_fn}) if $o{q_alu_fn};
	add_alu_info_to_trx_hash($self->{h_trx_hash}, $o{h_alu_fn}) if $o{h_alu_fn};
	
	push @{$self->{q_gene_hash}->{$_->{GENE_ID}}}, $_ foreach values %{$self->{q_trx_hash}};
	push @{$self->{h_gene_hash}->{$_->{GENE_ID}}}, $_ foreach values %{$self->{h_trx_hash}};
	
	return $self;
}

###
# PRIVATE METHODS

sub _get_blast_sites_for_trx
{
	my $self = shift;
	my($q_trx_id, $seed_len) = @_;
	
	my $q_trx = $self->{q_trx_hash}->{$q_trx_id};
	my $q_seq = $q_trx->{SEQ_ALU} || $q_trx->{SEQ};
	my $cmd   = get_blastn_command($self->{bdb_path},
		seed_len    => $seed_len,
		num_threads => $self->{num_threads},
		identity    => 70,
		outfmt      => 7,
	);
	
	my $pid = open2(*CHLD_OUT, *CHLD_IN, $cmd);
	print CHLD_IN ">$q_trx->{TRX_ID}\n$q_seq\n";
	close CHLD_IN;
	
	# Since we only searching the 'minus' strand the hit-start and hit-end are reverse in the BLASTn output -- let's fix it right away
#	my @header = qw(Q_TRX_ID H_TRX_ID IDENTITY ALI_LEN MISMATCHES GAPS_OPENINGS Q_START Q_END H_START H_END EVALUE SCORE);
	my @header = qw(Q_TRX_ID H_TRX_ID IDENTITY ALI_LEN MISMATCHES GAPS_OPENINGS Q_START Q_END H_END H_START EVALUE SCORE);
	my($all_hits, $n_lines) = ([], 0);
	while( <CHLD_OUT> )
	{
		$n_lines++;
		next if /^#/;
		
		my @vals = split /\s+/;
		warn "Wrong BLASTn line '$_'" and next if scalar(@vals) != scalar(@header);
		
		my %hit  = map { $header[$_] => $vals[$_] } 0..$#header;
		push @$all_hits, \%hit;
	}
	close CHLD_OUT;
	
	waitpid( $pid, 0 );   # to avoid zombie processes
	
	warn "The BLASTn output was empty for the query $q_trx_id" if $n_lines == 0;	
	
	return $all_hits;
}

sub _get_blast_site_gc_content
{
	my $self = shift;
	my($site) = @_;
	
	my $q_seq = $self->{q_trx_hash}->{$site->{Q_TRX_ID}}->{SEQ};
	my $h_seq = $self->{h_trx_hash}->{$site->{H_TRX_ID}}->{SEQ};
	die "Unknown TRX_ID=$site->{Q_TRX_ID} or $site->{H_TRX_ID}" if !$q_seq || !$h_seq;
	
	my $q_site_seq = substr($q_seq, $site->{Q_START}-1, $site->{Q_END}-$site->{Q_START}+1);
	my $h_site_seq = substr($h_seq, $site->{H_START}-1, $site->{H_END}-$site->{H_START}+1);
	my $site_seq   = uc($q_site_seq . $h_site_seq);
	my $num_gc     = $site_seq =~ tr/CG/CG/;
	
	return sprintf('%.1f', 100*$num_gc/length($site_seq));
	
}

sub _get_long_blast_sites_for_gene
{
	my $self = shift;
	my($q_gene_id, $seed_len, %opts) = @_;
	
	my $min_site_len    = $opts{min_site_len}    || undef;
	my $filter_overlaps = $opts{filter_overlaps} || 0;
	
	my $blast_sites = [];
	foreach my $q ( @{$self->{q_gene_hash}->{$q_gene_id}} )
	{
		print STDERR "[".localtime()."] BLASTing query gene $q_gene_id ($q->{TRX_ID})...            " if $self->{verbose};
		my $site_arr = $self->_get_blast_sites_for_trx($q->{TRX_ID}, $seed_len);
		
		# Add GENE_IDs and FILTER BY DUPLEX LENGTH!!!
		my $long_sites = [];
		foreach my $site ( sort {$b->{ALI_LEN} <=> $a->{ALI_LEN}} @$site_arr )
		{
			$site->{GC_CONTENT} = $self->_get_blast_site_gc_content($site);
			my $len_thr = $min_site_len || get_blast_site_len_thr($site->{GC_CONTENT}, $site->{IDENTITY});
			next if $site->{ALI_LEN} < $len_thr;
			
			next if $filter_overlaps && has_overlap_with_other_sites($site, $long_sites);
			
			$site->{ID}        = $self->{unique_id}++;
			$site->{Q_GENE_ID} = $self->{q_trx_hash}->{$site->{Q_TRX_ID}}->{GENE_ID};
			$site->{H_GENE_ID} = $self->{h_trx_hash}->{$site->{H_TRX_ID}}->{GENE_ID};
			$site->{Q_TRX_LEN} = length $self->{q_trx_hash}->{$site->{Q_TRX_ID}}->{SEQ};
			$site->{H_TRX_LEN} = length $self->{h_trx_hash}->{$site->{H_TRX_ID}}->{SEQ};
			push @$long_sites, $site;
		}
		
		my $num_sites = scalar(@$long_sites);
		my $num_genes = scalar unique(map { $_->{H_GENE_ID} } @$long_sites);
		
		print STDERR "$num_genes genes ($num_sites long sites) found\n" if $self->{verbose};
		push @$blast_sites, @$long_sites;
	}
	
	return $blast_sites;
}

###
# Each element of the $blast_sites has the following keys:
# Q_TRX_ID H_TRX_ID IDENTITY ALI_LEN MISMATCHES GAPS_OPENINGS Q_START Q_END H_START H_END EVALUE SCORE
# 
sub _filter_sites_by_isoforms
{
	my $self = shift;
	my($all_sites, $q_gene_id, $h_gene_id) = @_;
	
	my @gpair_sites = grep { $_->{Q_GENE_ID} eq $q_gene_id && $_->{H_GENE_ID} eq $h_gene_id } @$all_sites;
	
	# Grouping sites by the location on the query gives us HIT GROUPS
	my $site_to_h_trxs = {};
	foreach my $q_trx_id ( unique(map { $_->{Q_TRX_ID} } @gpair_sites) )
	{
		my $q_trx    = $self->{q_trx_hash}->{$q_trx_id};
		my @site_arr = grep { $_->{Q_TRX_ID} eq $q_trx->{TRX_ID} } @gpair_sites;
		foreach my $sgroup ( @{group_sites_by_subject(\@site_arr, length($q_trx->{SEQ}), 'Q')} )
		{
			my @h_trx_ids = unique( map { $_->{H_TRX_ID} } @$sgroup );
			foreach my $site ( @$sgroup )
			{
				warn "Site is present in more than one group: ".Dumper($site_to_h_trxs, $site) if $site_to_h_trxs->{$site->{ID}};
				$site_to_h_trxs->{ $site->{ID} } = \@h_trx_ids
			}
		}
	}
	
	# Grouping sites by the location on each hit gives us QUERY GROUPS
	my $site_to_q_trxs = {};
	foreach my $h_trx_id ( unique(map { $_->{H_TRX_ID} } @gpair_sites) )
	{
		my $h_trx    = $self->{h_trx_hash}->{$h_trx_id};
		my @site_arr = grep { $_->{H_TRX_ID} eq $h_trx->{TRX_ID} } @gpair_sites;
		foreach my $sgroup ( @{group_sites_by_subject(\@site_arr, length($h_trx->{SEQ}), 'H')} )
		{
			my @q_trx_ids = unique( map { $_->{Q_TRX_ID} } @$sgroup );
			foreach my $site ( @$sgroup )
			{
				warn "Site is present in more than one group: ".Dumper($site_to_h_trxs, $site) if $site_to_q_trxs->{$site->{ID}};
				$site_to_q_trxs->{ $site->{ID} } = \@q_trx_ids
			}
		}
	}
	
	my $good_sites = [];
	foreach my $site ( @gpair_sites )
	{
		my $num_q   = scalar @{$site_to_q_trxs->{$site->{ID}}};
		my $num_h   = scalar @{$site_to_h_trxs->{$site->{ID}}};
		my $total_q = scalar @{$self->{q_gene_hash}->{$q_gene_id}};
		my $total_h = scalar @{$self->{h_gene_hash}->{$h_gene_id}};
		push @$good_sites, $site if $num_q == $total_q && $num_h == $total_h;
	}
	
	return $good_sites;
}

sub _create_blast_site_files_with_flanks
{
	my $self = shift;
	my($flank_len, $dir) = @_;
	
	mkdir($dir) if !-d $dir;
	system("chmod 0777 $dir");
	
	my $blast_gpairs = $self->{blast_sites};
	my $site_arr     = [];
	foreach my $q_gene_id ( keys %$blast_gpairs )
	{
		foreach my $h_gene_id ( keys %{$blast_gpairs->{$q_gene_id}} )
		{
			foreach my $site ( values %{$blast_gpairs->{$q_gene_id}{$h_gene_id}} )
			{
				my %h = %$site;
				
				my $q_seq  = $self->{q_trx_hash}->{$site->{Q_TRX_ID}}->{SEQ};
				my $h_seq  = $self->{h_trx_hash}->{$site->{H_TRX_ID}}->{SEQ};
				my $q_site = _get_site_with_flanks($q_seq, $site->{Q_START}, $site->{Q_END}, $flank_len);
				my $h_site = _get_site_with_flanks($h_seq, $site->{H_START}, $site->{H_END}, $flank_len);
				
				$h{Q_SITE_START} = $q_site->{START};
				$h{Q_SITE_END}   = $q_site->{START} + length($q_site->{SEQ}) - 1;
				$h{Q_SITE_FNA}   = "$dir/site$site->{ID}.q_seq.fna";
				
				$h{H_SITE_START} = $h_site->{START};
				$h{H_SITE_END}   = $h_site->{START} + length($h_site->{SEQ}) - 1;
				$h{H_SITE_FNA}   = "$dir/site$site->{ID}.h_seq.fna";
				
				warn "Something is wrong with SITE_IDs" and next if -e $h{Q_SITE_FNA} || -e $h{H_SITE_FNA};
				write_to_file(">q\n$q_site->{SEQ}\n", $h{Q_SITE_FNA});
				write_to_file(">h\n$h_site->{SEQ}\n", $h{H_SITE_FNA});
				
				push @$site_arr, \%h;
			}
		}
	}
	
	return $site_arr;
}

sub _compute_coexpression_in_threads
{
	my $self = shift;
	my($dir, $bash_prefix) = @_;
	
	my $R_script = "$self->{tmp_dir}/coexpression.R";
	_create_coexpression_R( $R_script );
	
	my($gpairs, $txt_arr) = ($self->{blast_sites}, []);
	foreach my $q_gene_id ( sort keys %$gpairs )
	{
		my $q_fn = "$dir/$q_gene_id";
		if( !-e $q_fn )
		{
			warn "\nExpression profile file '$q_fn' doesn't exist!!!";
			$self->remove_gene_pair($q_gene_id, $_) foreach keys %{$gpairs->{$q_gene_id}};
			next;
		}
		
		foreach my $h_gene_id ( sort keys %{$gpairs->{$q_gene_id}} )
		{
			my $h_fn = "$dir/$h_gene_id";
			if( !-e $h_fn )
			{
				warn "\nExpression profile file '$h_fn' doesn't exist!!!";
				$self->remove_gene_pair($q_gene_id, $h_gene_id);
				next;
			}
			
			push @$txt_arr, "Rscript $R_script $dir $q_gene_id $h_gene_id";
			$txt_arr->[-1] = "echo 'Computing correlation for $q_gene_id - $h_gene_id...' 1>&2\n$txt_arr->[-1]" if $self->{verbose};
		}
	}
	
	my $size = int(scalar(@$txt_arr)/$self->{num_threads});
	$size++ if scalar(@$txt_arr) % $self->{num_threads};
	
	my($bash_arr, $out_files) = ([], []);
	foreach my $i (1..$self->{num_threads})
	{
		last if scalar(@$txt_arr) == 0;
		my @todo    = $i == $self->{num_threads} ? @$txt_arr : splice(@$txt_arr,0,$size);
		my $bash_fn = "$self->{tmp_dir}/$bash_prefix.part$i.sh";
		write_to_file( join("\n",@todo) , $bash_fn);
		chmod 0755, $bash_fn;
		
		my $out_fn = "$self->{tmp_dir}/$bash_prefix.part$i.txt";
		push @$bash_arr, "$bash_fn > $out_fn";
		push @$out_files, $out_fn;
	}
	
	run_commands_in_threads($bash_arr);
	
	return $out_files;
}

sub _run_bifold_in_threads
{
	my $self = shift;
	my($site_files, $dir, $bash_prefix, %opts) = @_;
	
	my $n_structs = $opts{num_structs} || $self->{num_structs};
	# http://stackoverflow.com/a/876242/310453
	my $redirect  = $opts{bifold_verbose} ? '1>&2' : '> /dev/null 2>&1';
	my $txt_arr   = [];
	foreach my $site ( @$site_files )
	{
		$site->{Q_FOLD_CT} = "$dir/site$site->{ID}.q_fold.ct";
		push @$txt_arr, "Fold $site->{Q_SITE_FNA} $site->{Q_FOLD_CT} --maximum 1 $redirect";
		
		$site->{H_FOLD_CT} = "$dir/site$site->{ID}.h_fold.ct";
		push @$txt_arr, "Fold $site->{H_SITE_FNA} $site->{H_FOLD_CT} --maximum 1 $redirect";
		
		$site->{BIFOLD_CT} = "$dir/site$site->{ID}.bifold.ct";
		push @$txt_arr, "bifold $site->{Q_SITE_FNA} $site->{H_SITE_FNA} $site->{BIFOLD_CT} --maximum $n_structs $redirect";
		
		$txt_arr->[-1] = "echo '[$bash_prefix] Running bifold for $site->{Q_GENE_ID}-$site->{H_GENE_ID}...' 1>&2\n$txt_arr->[-1]" if $self->{verbose};
	}
	
	my $size = int(scalar(@$txt_arr)/$self->{num_threads});
	$size++ if scalar(@$txt_arr) % $self->{num_threads};
	
	my $bash_arr = [];
	foreach my $i (1..$self->{num_threads})
	{
		last if scalar(@$txt_arr) == 0;
		my @todo    = $i == $self->{num_threads} ? @$txt_arr : splice(@$txt_arr,0,$size);
		my $bash_fn = "$self->{tmp_dir}/$bash_prefix.part$i.sh";
		write_to_file( join("\n",@todo) , $bash_fn);
		chmod 0755, $bash_fn;
		push @$bash_arr, $bash_fn;
	}
	
	run_commands_in_threads($bash_arr);
}

sub _read_bifold_files
{
	my $self = shift;
	my($site_files, $db_size) = @_;
	
	print STDERR "[".localtime()."] Reading bifold files...\n" if $self->{verbose};
	my $bifold_sites = {};
	foreach my $site ( @$site_files )
	{
		my $q_ct  = ASSA::Ct_formatManager->new( $site->{Q_FOLD_CT} );
		my $h_ct  = ASSA::Ct_formatManager->new( $site->{H_FOLD_CT} );
		my $bi_ct = ASSA::Ct_formatManager->new( $site->{BIFOLD_CT} );
		if( !$q_ct || !$h_ct || !$bi_ct )
		{
			warn "Couldn't read some of the .ct files: $site->{Q_FOLD_CT}, $site->{H_FOLD_CT}, $site->{BIFOLD_CT}!!!";
			$self->remove_blast_site($site->{Q_GENE_ID}, $site->{H_GENE_ID}, $site->{ID});
			next;
		}
		
		my $q_enrg     = $q_ct->get_full_energy;
		my $h_enrg     = $h_ct->get_full_energy;
		my $bi_enrg    = $bi_ct->get_full_energy;
		my $q_site_len = length $bi_ct->get_seq('query');
		my $h_site_len = length $bi_ct->get_seq('hit');
		my $dpl_len    = $bi_ct->get_total_dpl_ali_len();
		my $dpl_p      = compute_dpl_len_pvalue($self->{R}, $q_site_len, $h_site_len, $dpl_len);
		my $ddG        = sprintf('%.1f', $bi_enrg - ($q_enrg + $h_enrg));
		my $ddG_p      = compute_ddG_pvalue($self->{R}, $q_site_len, $h_site_len, $ddG);
		my $dpl_arr    = $bi_ct->dpl_arr(q_offset => $site->{Q_SITE_START}-1, h_offset => $site->{H_SITE_START}-1, uc_keys => 1);
		my $dpl_txt    = $bi_ct->get_dpl_txt(
			q_gene_id => $site->{Q_GENE_ID}, q_trx_id => $site->{Q_TRX_ID}, q_offset => $site->{Q_SITE_START}-1,
			h_gene_id => $site->{H_GENE_ID}, h_trx_id => $site->{H_TRX_ID}, h_offset => $site->{H_SITE_START}-1,
		);
		
		my $id = $self->{unique_id}++;
		$bifold_sites->{$site->{Q_GENE_ID}}->{$site->{H_GENE_ID}}->{$id} = {
			ID             => $id,                           DPL_LEN        => $dpl_len,
			BLAST_SITE_ID  => $site->{ID},                   DPL_LEN_PVALUE => $dpl_p,
			Q_TRX_ID       => $site->{Q_TRX_ID},             DDG            => $ddG,
			H_TRX_ID       => $site->{H_TRX_ID},             DDG_PVALUE     => $ddG_p,
			Q_SEQ          => $bi_ct->get_seq('query'),      
			Q_START        => $site->{Q_SITE_START},         DPL_TXT        => $dpl_txt,
			Q_END          => $site->{Q_SITE_END},           DPL_ARR        => $dpl_arr,
			H_SEQ          => $bi_ct->get_seq('hit'),
			H_START        => $site->{H_SITE_START},
			H_END          => $site->{H_SITE_END},
		};
	}
	
	$self->{bifold_sites} = $bifold_sites;
}

sub _blast_and_bifold_must_overlap
{
	my $self = shift;
	
	my $bifold_sites = $self->{bifold_sites};
	foreach my $q_gene_id ( keys %$bifold_sites )
	{
		foreach my $h_gene_id ( keys %{$bifold_sites->{$q_gene_id}} )
		{
			foreach my $bi_site ( values %{$bifold_sites->{$q_gene_id}{$h_gene_id}} )
			{
				# the original BLAST duplex must be reconstructed at least in part
				my $blast_site = $self->{blast_sites}{$q_gene_id}{$h_gene_id}{$bi_site->{BLAST_SITE_ID}};
				my $has_blast  = has_double_overlap_with_other_sites($blast_site, $bi_site->{DPL_ARR});
				$self->remove_bifold_site($q_gene_id, $h_gene_id, $bi_site->{ID}) if !$has_blast;
			}
		}
	}
}

sub _make_gene_pairs
{
	my $self = shift;
	my($evalue_thr, $db_size) = @_;
	
	my $bifold_sites = $self->{bifold_sites};
	$self->{gene_pairs} = {};
	foreach my $q_gene_id ( keys %$bifold_sites )
	{
		foreach my $h_gene_id ( keys %{$bifold_sites->{$q_gene_id}} )
		{
			my @all_sites = values %{$bifold_sites->{$q_gene_id}{$h_gene_id}};
			my $q_len     = length $self->{q_trx_hash}->{$all_sites[0]{Q_TRX_ID}}->{SEQ};
			my $h_len     = length $self->{h_trx_hash}->{$all_sites[0]{H_TRX_ID}}->{SEQ};
			my($score, $evalue);
			if( grep { $_->{DDG_PVALUE} == 0 } @all_sites )
			{
				($score, $evalue) = (999999,0);
			}
			else
			{
				$score  = sum( map { -1*log10($_->{DDG_PVALUE}) } @all_sites );
				$evalue = compute_ddG_evalue($q_len, $h_len, $score, $db_size);
			}
			
			if( $evalue > $evalue_thr )
			{
				$self->remove_gene_pair($q_gene_id,$h_gene_id);
				next;
			}
			
			$self->{gene_pairs}{$q_gene_id}{$h_gene_id} = {
				Q_TRX_ID   => $all_sites[0]{Q_TRX_ID},
				H_TRX_ID   => $all_sites[0]{H_TRX_ID},
				NUM_SITES  => scalar( @all_sites ),
				SCORE      => sprintf('%.2f', $score),
				DDG_EVALUE => $evalue,
			};
		}
	}
}

sub _remove_overlapping_bifold_sites
{
	my $self = shift;
	
	my $bifold_sites = $self->{bifold_sites};
	foreach my $q_gene_id ( keys %$bifold_sites )
	{
		foreach my $h_gene_id ( keys %{$bifold_sites->{$q_gene_id}} )
		{
			my @all_sites  = values %{$bifold_sites->{$q_gene_id}{$h_gene_id}};
			
			my $no_overlap = [];
			foreach my $site ( sort { $a->{DDG_PVALUE} <=> $b->{DDG_PVALUE} } @all_sites )
			{
				if( has_overlap_with_other_sites($site, $no_overlap) )
				{
					delete $bifold_sites->{$q_gene_id}->{$h_gene_id}->{$site->{ID}};
				}
				else
				{
					push @$no_overlap, $site;
				}
			}
		}
	}
}

###
# PUBLIC METHODS

sub has_coexpression{ return $_[0]->{coexpression} ? 1 : 0; }
sub has_bifold      { return $_[0]->{bifold_sites} ? 1 : 0; }

sub run_bifold
{
	my $self = shift;
	my($flank_len, $name, %opts) = @_;
	
	my $pvalue_ddg      = defined $opts{pvalue_ddg}      ? $opts{pvalue_ddg}      : 1.0;
	my $evalue_ddg      = defined $opts{evalue_ddg}      ? $opts{evalue_ddg}      : 1.0;
	my $db_size         = defined $opts{db_size}         ? $opts{db_size}         : scalar(keys %{$self->{h_gene_hash}});
	my $filter_overlaps = defined $opts{filter_overlaps} ? $opts{filter_overlaps} : 1;    # there must not be two overlapping bifold regions
	my $filter_bbo      = defined $opts{filter_bbo}      ? $opts{filter_bbo}      : 1;    # at least one bifold duplex must overlap with BLAST site
	
	$self->_merge_overlapping_blast_sites_with_flank( $flank_len );
	
	my $short_dir   = "$self->{tmp_dir}/$name";
	my $short_sites = $self->_create_blast_site_files_with_flanks($flank_len, $short_dir);
	$self->_run_bifold_in_threads($short_sites, $short_dir, $name, %opts);
	$self->_read_bifold_files($short_sites, $db_size);                                    # create $self->{bifold_sites}
	$self->_blast_and_bifold_must_overlap()   if $filter_bbo;
	$self->_remove_overlapping_bifold_sites() if $filter_overlaps;
	$self->_make_gene_pairs($evalue_ddg, $db_size);
	
	if( $self->{out_dir} )
	{
		my $txt = '';
		foreach my $q_gene_id ( sort keys %{$self->{gene_pairs}} )
		{
			foreach my $h_gene_id ( keys %{$self->{gene_pairs}->{$q_gene_id}} )
			{
				my $gpair = $self->{gene_pairs}{$q_gene_id}{$h_gene_id};
				$txt .= "$q_gene_id\t$h_gene_id\t$gpair->{SCORE}\t$gpair->{DDG_EVALUE}\n";
			}
		}
		write_to_file($txt, "$self->{out_dir}/hits_after_bifold_$flank_len.txt");
	}
	
	print STDERR "[".localtime()."] Number of gene pairs after bifold filter ($name) = ".$self->get_num_gene_pairs."\n" if $self->{verbose};
}

sub run_blast
{
	my $self = shift;
	my($seed_len, %opts) = @_;
	
	# makeblastdb
	open(PIPE, "| makeblastdb -dbtype nucl -out $self->{bdb_path} -title 'tmp_db'> /dev/null 2>&1") or die "Can't open BLAST pipe: $!";
	print PIPE ">$_->{TRX_ID}\n$_->{SEQ}\n" foreach values %{$self->{h_trx_hash}};
	close PIPE;
	
	my($gene_pair_hash, $total_gpairs, $total_sites) = ({}, 0, 0);
	foreach my $q_gene_id ( sort keys %{$self->{q_gene_hash}} )
	{
		# Apply MIN_BLAST_SITE_LEN
		my $all_sites = $self->_get_long_blast_sites_for_gene($q_gene_id, $seed_len, %opts);
		
		# Apply FILTER: 100% interacting isoforms
		$gene_pair_hash->{$q_gene_id} = {};

# Uncomment this and comment out the foreach look to speed up the process
#$gene_pair_hash->{$q_gene_id}->{$_->{H_GENE_ID}}->{$_->{ID}} = $_ foreach @$all_sites;
#$total_gpairs = scalar( unique( map { $_->{H_GENE_ID} } @$all_sites ) );
#$total_sites  = scalar( @$all_sites );

		foreach my $h_gene_id ( unique( map { $_->{H_GENE_ID} } @$all_sites ) )
		{
			my $gpair_sites = [grep { $_->{Q_GENE_ID} eq $q_gene_id && $_->{H_GENE_ID} eq $h_gene_id } @$all_sites];
			   $gpair_sites = $self->_filter_sites_by_isoforms($gpair_sites, $q_gene_id, $h_gene_id) if $opts{filter_isoforms};
			my $best_sites  = get_best_trx_sites( $gpair_sites );
			
			$gene_pair_hash->{$q_gene_id}->{$h_gene_id} = $best_sites if %$best_sites;
			
			$total_gpairs++;
			$total_sites += scalar(keys %$best_sites);
		}
	}
	
	if( $self->{out_dir} )
	{
		my $txt = '';
		foreach my $q_gene_id ( sort keys %$gene_pair_hash )
		{
			$txt .= $_ foreach map { "$q_gene_id\t$_\n" } keys %{$gene_pair_hash->{$q_gene_id}};
		}
		write_to_file($txt, "$self->{out_dir}/hits_after_blast.txt");
	}
	
	$self->{blast_sites} = $gene_pair_hash;
	$self->_merge_overlapping_blast_sites_with_flank( 0 );
	
	print STDERR "[".localtime()."] Number of gene pairs after BLAST search = $total_gpairs ($total_sites sites)\n"    if $self->{verbose};
	print STDERR "[".localtime()."] Number of gene pairs after 100% isoforms filter = ".$self->get_num_gene_pairs."\n" if $self->{verbose} && $opts{filter_isoforms};
}

sub remove_blast_site
{
	my $self = shift;
	my($q_gene_id, $h_gene_id, $blast_id) = @_;
	
	delete $self->{blast_sites}{$q_gene_id}{$h_gene_id}{$blast_id};
	
	$self->remove_gene_pair($q_gene_id,$h_gene_id) if !keys( %{$self->{blast_sites}{$q_gene_id}{$h_gene_id}} );
}

sub remove_bifold_site
{
	my $self = shift;
	my($q_gene_id, $h_gene_id, $bifold_id) = @_;
	
	my $bifold_site = $self->{bifold_sites}{$q_gene_id}{$h_gene_id}{$bifold_id};
	warn "Unknown bifold site id '$bifold_id'" and return if !$bifold_site;
	
	$self->remove_blast_site($q_gene_id, $h_gene_id, $bifold_site->{BLAST_SITE_ID});
	delete $self->{bifold_sites}{$q_gene_id}{$h_gene_id}{$bifold_id};
	
	$self->remove_gene_pair($q_gene_id,$h_gene_id) if !keys( %{$self->{bifold_sites}{$q_gene_id}{$h_gene_id}} );
}

sub remove_gene_pair
{
	my $self = shift;
	my($q_gene_id,$h_gene_id) = @_;
	
	if( my $gpairs = $self->{blast_sites} )
	{
		delete $gpairs->{$q_gene_id}{$h_gene_id};
		delete $gpairs->{$q_gene_id} if !keys( %{$gpairs->{$q_gene_id}} );   # remove the empty q_gene_id key if needed
	}
	
	if( my $gpairs = $self->{bifold_sites} )
	{
		delete $gpairs->{$q_gene_id}{$h_gene_id};
		delete $gpairs->{$q_gene_id} if !keys( %{$gpairs->{$q_gene_id}} );   # remove the empty q_gene_id key if needed
	}
	
	if( my $gpairs = $self->{coexpression} )
	{
		delete $gpairs->{$q_gene_id}{$h_gene_id};
		delete $gpairs->{$q_gene_id} if !keys( %{$gpairs->{$q_gene_id}} );   # remove the empty q_gene_id key if needed
	}
}

sub remove_gene_pairs_with_alu
{
	my $self = shift;
	
	my $gene_pairs = $self->{blast_sites};
	foreach my $q_gene_id ( keys %$gene_pairs )
	{
		my $q_has_alu = grep { $_->{NUM_ALU} } @{$self->{q_gene_hash}->{$q_gene_id}};
		next if !$q_has_alu;
		
		foreach my $h_gene_id ( keys %{$gene_pairs->{$q_gene_id}} )
		{
			my $h_has_alu = grep { $_->{NUM_ALU} } @{$self->{h_gene_hash}->{$h_gene_id}};
			$self->remove_gene_pair($q_gene_id,$h_gene_id) if $q_has_alu && $h_has_alu;
		}
	}
	print STDERR "[".localtime()."] Number of gene pairs after Alu filter = ".$self->get_num_gene_pairs."\n" if $self->{verbose};
}

sub remove_cis_gene_pairs
{
	my $self = shift;
	
	my $gene_pairs = $self->{blast_sites};
	foreach my $q_gene_id ( keys %$gene_pairs )
	{
		foreach my $q_trx ( @{$self->{q_gene_hash}->{$q_gene_id}} )
		{
			my $q_trx_id  = $q_trx->{TRX_ID};
			my $cis_sites = $self->_get_blast_sites_for_trx($q_trx_id, $self->{cis_seed});
			foreach my $site ( @$cis_sites )
			{
				my $h_gene_id = $self->{h_trx_hash}->{$site->{H_TRX_ID}}->{GENE_ID};
				$self->{q_trx_hash}{$q_trx_id}{CIS_GENES}{$h_gene_id} = 1;
				$self->remove_gene_pair($q_gene_id,$h_gene_id);
			}
		}
	}
	print STDERR "[".localtime()."] Number of gene pairs after cis-filter = ".$self->get_num_gene_pairs."\n" if $self->{verbose};
}

sub filter_gene_pairs_by_coexpression
{
	my $self = shift;
	my($dir, $pvalue_pcc, %opts) = @_;
	
	my $out_files = $self->_compute_coexpression_in_threads(abs_path($dir), 'run_coexpression');
	
	$self->{coexpression} = {};
	foreach my $out_fn ( @$out_files )
	{
		foreach my $arr ( @{read_table_wo_header($out_fn)} )
		{
			my($q_gene_id, $h_gene_id, $num, $pcc, $pvalue) = @$arr;
			if( $pvalue < $pvalue_pcc )
			{
				$self->{coexpression}->{$q_gene_id}->{$h_gene_id} = {
					PCC        => sprintf('%.3f', $pcc   ),
					PCC_PVALUE => sprintf('%.2e', $pvalue),
					PCC_NUM    => $num,
				};
			}
			else
			{
				$self->remove_gene_pair($q_gene_id,$h_gene_id);
			}
		}
	}
	
	print STDERR "[".localtime()."] Number of gene pairs after coexpression filter = ".$self->get_num_gene_pairs."\n" if $self->{verbose};
}

sub get_num_gene_pairs
{
	my $self = shift;
	
	my($gpair_hash, $num_gpairs) = ($self->{blast_sites}, 0);
	foreach my $q_gene_id ( keys %$gpair_hash )
	{
		$num_gpairs += scalar keys %{$gpair_hash->{$q_gene_id}};
	}
	
	return $num_gpairs;
}

sub get_notes_for_q_gene
{
	my $self = shift;
	my($q_gene_id) = @_;
	
	my $q_notes = [];
	push @$q_notes, 'alu' if grep { $_->{NUM_ALU}   } @{$self->{q_gene_hash}{$q_gene_id}};
	push @$q_notes, 'cis' if grep { $_->{CIS_GENES} } @{$self->{q_gene_hash}{$q_gene_id}};
	
	return join(';', @$q_notes),
}

sub get_blast_site_by_id
{
	my $self = shift;
	my($id) = @_;
	
	my $gpair_hash = $self->{blast_sites};
	foreach my $q_gene_id ( keys %$gpair_hash )
	{
		foreach my $h_gene_id ( keys %{$gpair_hash->{$q_gene_id}} )
		{
			foreach my $site ( values %{$gpair_hash->{$q_gene_id}{$h_gene_id}} )
			{
				return $site if $site->{ID} eq $id;
			}
		}
	}
	
	return undef;
}

sub get_blast_site_aoh
{
	my $self = shift;
	
	my($site_hash, $site_aoh) = ($self->{blast_sites}, []);
	foreach my $q_gene_id ( keys %$site_hash )
	{
		foreach my $h_gene_id ( keys %{$site_hash->{$q_gene_id}} )
		{
			push @$site_aoh, values %{$site_hash->{$q_gene_id}{$h_gene_id}};
		}
	}
	
	return $site_aoh;
}

sub get_bifold_site_aoh
{
	my $self = shift;
	
	my($site_hash, $site_aoh) = ($self->{bifold_sites}, []);
	foreach my $q_gene_id ( keys %$site_hash )
	{
		foreach my $h_gene_id ( keys %{$site_hash->{$q_gene_id}} )
		{
			my $coexpression = $self->{coexpression} ? $self->{coexpression}{$q_gene_id}{$h_gene_id} : {};
			foreach my $site ( values %{$site_hash->{$q_gene_id}{$h_gene_id}} )
			{
				my %hash = %$site;
				$hash{COEXPRESSION} = $coexpression;
				$hash{Q}            = $self->{q_trx_hash}->{ $site->{Q_TRX_ID} };
				$hash{H}            = $self->{h_trx_hash}->{ $site->{H_TRX_ID} };
				push @$site_aoh, \%hash;
			}
		}
	}
	
	return [ sort { $a->{DDG_PVALUE} <=> $b->{DDG_PVALUE} } @$site_aoh ];
}

sub get_gene_pair_aoh
{
	my $self = shift;
	
	my($gpair_hash, $gpair_aoh) = ($self->{bifold_sites}, []);
	foreach my $q_gene_id ( keys %$gpair_hash )
	{
		foreach my $h_gene_id ( keys %{$gpair_hash->{$q_gene_id}} )
		{
			my @all_sites = sort { $a->{DDG_PVALUE} <=> $b->{DDG_PVALUE} } values %{$gpair_hash->{$q_gene_id}{$h_gene_id}};
#			my %best_site = %{$all_sites[0]};
			
			my %best_site = ();
			my $q_trx_id  = $self->{gene_pairs}{$q_gene_id}{$h_gene_id}{Q_TRX_ID};
			my $h_trx_id  = $self->{gene_pairs}{$q_gene_id}{$h_gene_id}{H_TRX_ID};
			$best_site{Q}            = $self->{q_trx_hash}->{$q_trx_id};
			$best_site{H}            = $self->{h_trx_hash}->{$h_trx_id};
			$best_site{COEXPRESSION} = $self->{coexpression} ? $self->{coexpression}{$q_gene_id}{$h_gene_id} : {};
			$best_site{BIFOLD_SITES} = \@all_sites;
			$best_site{NUM_SITES}    = $self->{gene_pairs}{$q_gene_id}{$h_gene_id}{NUM_SITES};
			$best_site{DPL_LEN}      = sum( map { $_->{DPL_LEN} } @all_sites );
			$best_site{SCORE}        = $self->{gene_pairs}{$q_gene_id}{$h_gene_id}{SCORE};
			$best_site{DDG_EVALUE}   = $self->{gene_pairs}{$q_gene_id}{$h_gene_id}{DDG_EVALUE};
			
			push @$gpair_aoh, \%best_site;
		}
	}
	
	return [ sort { $a->{DDG_EVALUE} <=> $b->{DDG_EVALUE} } @$gpair_aoh ];
}

###
# SUBROUTINES

sub run_commands_in_threads
{
	my($bash_arr) = @_;
	
	# http://stackoverflow.com/a/9253948
	my $all_pids = {};
	foreach my $bash_fn ( @$bash_arr )
	{
		my $pid = fork();
		$all_pids->{$pid} = $bash_fn;
		
		if( $pid == 0 )    # child
		{
			exec( $bash_fn );
			die "Exec $bash_fn failed: $!\n";
		}
		elsif( !defined $pid )
		{
			warn "Fork $bash_fn failed: $!\n";
		}
	}
	
	# http://stackoverflow.com/a/16980565
	while( keys %$all_pids )
	{
		#warn Dumper($all_pids);
		foreach my $pid ( keys %$all_pids )
		{
			if( waitpid($pid, WNOHANG) != 0 )
			{
				# Child is done
				delete $all_pids->{$pid};
			}
		}
		sleep(1);
	}
}

sub _create_coexpression_R
{
	my($fn) = @_;
	
	write_to_file(q~
# http://stackoverflow.com/a/2154190
args     <- commandArgs(trailingOnly = TRUE)
dir      <- args[1]
q.name   <- args[2]
h.name   <- args[3]
min.expr <- 1

# Useful graph: http://stats.stackexchange.com/questions/17371/example-of-strong-correlation-coefficient-with-a-high-p-value
min.vals <- 5

q.fn   <- paste(dir, q.name, sep = "/")
h.fn   <- paste(dir, h.name, sep = "/")
q.vals <- read.table(q.fn, header=FALSE)
h.vals <- read.table(h.fn, header=FALSE)

# make data.frame and filter numbers by min_expression
df <- data.frame(q.vals, h.vals)
colnames(df) <- c("q.vals", "h.vals")
look.at <- df$q.vals > min.expr & df$h.vals > min.expr

if( sum(look.at) > min.vals ) {
	pcc <- cor.test(df$q.vals[look.at] , df$h.vals[look.at] , method = "pearson")
	cat(q.name, h.name, sum(look.at), pcc$estimate, pcc$p.value, "\n", sep = "\t")
} else {
	cat(q.name, h.name, sum(look.at), 0,            0,           "\n", sep = "\t")
}
~, $fn);
}

sub _get_site_with_flanks
{
	my($seq, $start, $end, $flank_len) = @_;
	
	my $offset = $start - $flank_len - 1;
	$offset = 0 if $offset < 0;   # start cutting at the start of the sequence
	
	my $site_len = 2*$flank_len + ($end - $start + 1);
	my $site_seq = '';
	if( $offset + $site_len > length($seq) )
	{
		# stop cutting at the end of the sequence
		$site_seq = substr $seq, $offset;
	}
	else
	{
		$site_seq = substr $seq, $offset, $site_len;
	}
	
	return {START => $offset+1, SEQ => uc($site_seq)};
}

# GC = [0..100], cmpl = [0..100]
sub get_blast_site_len_thr
{
	my($gc, $cmpl) = @_;
	return 462.276074461 - 0.307494149*$gc - 8.543762823*$cmpl + 0.001881507*$gc**2 + 0.041391051*$cmpl**2;
}

sub compute_ddG_pvalue
{
	my($R, $len1, $len2, $ddG) = @_;
	
	my $short = $len1 < $len2 ? $len1 : $len2;
	my $long  = $len1 < $len2 ? $len2 : $len1;
	my $mean  = -9.1129539 - 0.0041405*$short - 0.0021826*$long;
	my $sd    =  4.3317524 + 0.0006908*$short + 0.0002260*$long;
	
	return sprintf '%.2e', ASSA::Gauss::cdf($ddG, $mean, $sd);
}

sub compute_dpl_len_pvalue
{
	my($R, $len1, $len2, $dpl_len) = @_;
	
	my $short = $len1 < $len2 ? $len1 : $len2;
	my $long  = $len1 < $len2 ? $len2 : $len1;
	my $mean  = -1.586727 + 0.070465*$short + 0.036295*$long;
	my $sd    = 13.933122 + 0.031977*$short + 0.014195*$long;
	
	return sprintf('%.2e', 1 - ASSA::Gauss::cdf($dpl_len, $mean, $sd));
}

sub compute_ddG_evalue
{
	my($len1, $len2, $score, $db_size) = @_;
	
	my $m      = $len1 / 1000;        # length should be in kb
	my $n      = $len2 / 1000;        # length should be in kb
	my $K      = 0.11385;
	my $lambda = 1.166;
	
	return sprintf '%.2e', $db_size * $K * $m * $n * exp(-1 * $lambda * $score);
}

sub has_overlap_with_other_sites
{
	my($site, $site_arr) = @_;
	
	# check overlaps on the query OR on the hit
	foreach my $s ( @$site_arr )
	{
		next if $s->{Q_TRX_ID} ne $site->{Q_TRX_ID} || $s->{H_TRX_ID} ne $site->{H_TRX_ID};
		return 1 if intersect_digits($site->{Q_START},$site->{Q_END},$s->{Q_START},$s->{Q_END});
		return 1 if intersect_digits($site->{H_START},$site->{H_END},$s->{H_START},$s->{H_END});
	}
	
	return 0;
}

sub has_double_overlap_with_other_sites
{
	my($site, $site_arr) = @_;
	
	# check overlaps on the query AND on the hit
	foreach my $s ( @$site_arr )
	{
		if( $s->{Q_TRX_ID} && $site->{Q_TRX_ID} && $s->{H_TRX_ID} && $site->{H_TRX_ID} )
		{
			next if $s->{Q_TRX_ID} ne $site->{Q_TRX_ID} || $s->{H_TRX_ID} ne $site->{H_TRX_ID};
		}
		
		if( intersect_digits($site->{Q_START},$site->{Q_END},$s->{Q_START},$s->{Q_END})
		    &&
		    intersect_digits($site->{H_START},$site->{H_END},$s->{H_START},$s->{H_END}))
		{
			return $s;
		}
	}
	
	return undef;
}

sub _merge_overlapping_blast_sites_with_flank
{
	my $self = shift;
	my($flank_len) = @_;
	
	# group sites by trx pair
	my($trx_pairs, $all_sites) = ({}, $self->get_blast_site_aoh);
	push @{$trx_pairs->{"$_->{Q_TRX_ID}::$_->{H_TRX_ID}"}}, $_ foreach @$all_sites;
	
	my $merged_sites = [];
	foreach my $site_arr ( values %$trx_pairs )
	{
		my $q_trx_id = $site_arr->[0]{Q_TRX_ID};
		my $h_trx_id = $site_arr->[0]{H_TRX_ID};
		my $q_len    = length($self->{q_trx_hash}{$q_trx_id}{SEQ});
		my $h_len    = length($self->{h_trx_hash}{$h_trx_id}{SEQ});
		
		my $check_for_overlaps = 1;
		while( $check_for_overlaps )
		{
			$check_for_overlaps = 0;
			
			my $no_overlap = [];
			foreach my $site ( sort { $a->{Q_START} <=> $b->{Q_START} } @$site_arr )
			{
				my %site_with_flanks = %$site;
				$site_with_flanks{Q_START} = max(1,      $site->{Q_START} - 2*$flank_len);
				$site_with_flanks{Q_END}   = min($q_len, $site->{Q_END}   + 2*$flank_len);
				$site_with_flanks{H_START} = max(1,      $site->{H_START} - 2*$flank_len);
				$site_with_flanks{H_END}   = min($h_len, $site->{H_END}   + 2*$flank_len);
				if( my $ovlp_site = has_double_overlap_with_other_sites(\%site_with_flanks, $no_overlap) )
				{
					$check_for_overlaps = 1;
					
					# Extend the found overlapping site
					$ovlp_site->{Q_START} = min($ovlp_site->{Q_START}, $site->{Q_START});
					$ovlp_site->{Q_END}   = max($ovlp_site->{Q_END},   $site->{Q_END}  );
					$ovlp_site->{H_START} = min($ovlp_site->{H_START}, $site->{H_START});
					$ovlp_site->{H_END}   = max($ovlp_site->{H_END},   $site->{H_END}  );
				}
				else
				{
					push @$no_overlap, $site;
				}
			}
			$site_arr = $no_overlap;
		}
		push @$merged_sites, @$site_arr;
	}
	
	# save the merged sites to $self->{blast_sites}
	$self->{blast_sites} = {};
	$self->{blast_sites}{$_->{Q_GENE_ID}}{$_->{H_GENE_ID}}{$_->{ID}} = $_ foreach @$merged_sites;
}

sub get_best_trx_sites
{
	my($site_arr) = @_;
	
	# Determine the TRX pair that has the best site: best e-value, then min total len
	my($best) = sort {
	                    $b->{ALI_LEN}                    <=>  $a->{ALI_LEN}                      or
	                    $a->{EVALUE}                     <=>  $b->{EVALUE}                       or
	                   ($a->{Q_TRX_LEN}+$a->{H_TRX_LEN}) <=> ($b->{Q_TRX_LEN}+$b->{H_TRX_LEN})
	                 } @$site_arr;
	
	my $best_trx_sites = {};
	foreach my $site ( @$site_arr )
	{
		if( $site->{Q_TRX_ID} eq $best->{Q_TRX_ID} && $site->{H_TRX_ID} eq $best->{H_TRX_ID} )
		{
			warn "Something is wrong with the SITE_ID: ".Dumper($site) if $best_trx_sites->{$site->{ID}};
			$best_trx_sites->{$site->{ID}} = $site;
		}
	}
	
	return $best_trx_sites;
}

sub group_sites_by_subject
{
	my($site_arr, $subject_len, $prefix) = @_;
	
	my $subject_start = $prefix.'_START';
	my $subject_end   = $prefix.'_END';
	my $subject_str   = '.'x$subject_len;
	foreach my $site ( @$site_arr )
	{
		warn "Something is wrong for BLAST site: ".Dumper($site) and next if $site->{$subject_end} > $subject_len;
		
		my $offset  = $site->{$subject_start}-1;
		my $dpl_len = $site->{$subject_end} - $site->{$subject_start} + 1;
		my $dpl_str = 'D'x$dpl_len;
		
		my $prefix = substr($subject_str, 0, $offset);
		my $suffix = substr($subject_str, $offset+$dpl_len);
		$subject_str = $prefix . $dpl_str . $suffix;
	}
	die "The subject string is wrong: $subject_str" if length($subject_str) != $subject_len;
	
	# $subject_str = '...DDDDD.....DD...DDDD.....'   =>   numbers
	my $site_groups = [];
	while( $subject_str =~ /D+/g )
	{
		my $start = length($`)+1;
		my $len   = length($&);
		push @$site_groups, {start => $start, len => $len, end => $start+$len-1, sites => []};
	}
	
	foreach my $site ( @$site_arr )
	{
		my $group_found = 0;
		foreach my $sgroup ( @$site_groups )
		{
			if( $sgroup->{start} <= $site->{$subject_start} && $site->{$subject_end} <= $sgroup->{end} )
			{
				push @{$sgroup->{sites}}, $site;
				$group_found = 1;
				last;
			}
		}
		warn "Could not find group for site: ".Dumper($site_groups, $site) if !$group_found && @$site_groups;
	}
	
	return ah2a('sites', $site_groups);
}

sub add_alu_info_to_trx_hash
{
	my($trx_hash, $masked_fn) = @_;
	
	my $masked_hash = fasta2hash( $masked_fn );
	foreach my $trx_id ( keys %$masked_hash )
	{
		warn "Alu-masked file '$masked_fn' contains unknown trx '$trx_id'" and next if !$trx_hash->{$trx_id};
		
		my $seq_alu = $masked_hash->{$trx_id}{SEQ};
		my $num_alu = $seq_alu =~ s/(N+)/$1/g;
		$trx_hash->{$trx_id}{NUM_ALU} = $num_alu || 0;
		$trx_hash->{$trx_id}{SEQ_ALU} = $seq_alu if $num_alu;
	}
}

sub get_blastn_command
{
	my($db_path, %o) = @_;
	
	$o{outfmt}     = 6         if !defined $o{outfmt};
	$o{num_threads}= 8         if !defined $o{num_threads};
	$o{strand}     = 'minus'   if !        $o{strand};
	$o{n_hits}     = 100000000 if !        $o{n_hits};
	$o{evalue}     = 100000000 if !        $o{evalue};
	$o{seed_len}   = 12        if !        $o{seed_len};
	$o{identity}   = 70        if !        $o{identity};
	$o{allow_mism} = 1         if !defined $o{allow_mism};
	$o{allow_gaps} = 1         if !defined $o{allow_gaps};
	
	my $cmd = "blastn -task blastn -num_threads $o{num_threads} -outfmt $o{outfmt} -strand $o{strand}";
#	$cmd .=  ' -dust no -soft_masking false';
	$cmd .=  $db_path       ? " -db              $db_path"     : '';
	$cmd .=  $o{evalue}     ? " -evalue          $o{evalue}"   : '';
	$cmd .=  $o{seed_len}   ? " -word_size       $o{seed_len}" : '';
	$cmd .=  $o{identity}   ? " -perc_identity   $o{identity}" : '';
	$cmd .=  $o{query}      ? " -query           $o{query}"    : '';
	$cmd .=  $o{out}        ? " -out             $o{out}"      : '';
	$cmd .= !$o{allow_mism} ? ' -penalty         -1000'        : '';
	$cmd .= !$o{allow_gaps} ? ' -ungapped'                     : '';
	
	if( $o{outfmt} <= 3 )
	{
		# Warning: The parameter -max_target_seqs is ignored for output formats, 0,1,2,3.
		# Use -num_descriptions and -num_alignments to control output
		$cmd .= " -num_descriptions $o{n_hits} -num_alignments $o{n_hits}";
	}
	else
	{
		$cmd .= " -max_target_seqs $o{n_hits}";
	}
	
	return $cmd;
}

1;

