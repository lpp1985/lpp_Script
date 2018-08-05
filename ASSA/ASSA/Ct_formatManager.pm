package ASSA::Ct_formatManager;

use strict;
use warnings;

# $Id: Ct_formatManager.pm 2554 2016-04-21 08:07:35Z antonov $

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
#
# $self
#   |
#   |--->{max_structs}  --  [number] option
#   |--->{n_structs}    --  [number] number of structures that were read from the file
#   |--->{enrg}         --  [number] MFE energy
#   |
#   |--->{q_aoh}
#   |--->{h_aoh}
#   |     h_aoh[i]
#   |       |
#   |       |--->{i}
#   |       |--->{seq_i}
#   |       |--->{let}
#   |       |--->{paired}
#   |       |
#   |       |--->{unpaired}    --  [bool]
#   |       |--->{unpaired_p}  --  [0..1]
#   |       |
#   |       |--->{duplex}    --  [bool] inter-molecular duplex
#   |       |--->{duplex_i}  --  index of the h_aoh element that interacts with the nt in the q_aoh
#   |       |--->{duplex_p}  --  [probability: 0..1]
#   |       |
#   |       |--->{state}     --  secondary structure state: U (unpaired), D (duplex), S (secondary structure)
#   |
#   |--->{dpl_arr}
#   |     dpl_arr[i]
#   |        |
#   |        |--->{q_start_i}
#   |        |--->{q_end_i}
#   |        |--->{h_start_i}
#   |        |--->{h_end_i}
#   |        |--->{q_ali}
#   |        |--->{h_ali_r}
#   |        |--->{match_str}
#   |        |--->{n_matches}
#   |        |--->{identity}    --  [0..100]
#   |        |--->{len}         --  [bp]
#   |
#

use Data::Dumper;
use Carp qw(confess); 

use ASSA::Lib qw(revcomp read_table_wo_header_from_fh);

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw();

###
# CONSTANTS

# For the SW-algorithm
my ($S, $D, $L, $U) = qw(* \ - |);   # Start/Stop, Diagonal, Left, Up -- is used in @from_array

###
# CONSTRUCTOR
sub new
{
	my $class = shift;
	my($fn, %opts) = @_;
	
	if( !-e $fn )
	{
		warn "The file '$fn' does not exist!!!";
		return undef;
	}
	
	my $self = bless {
		enrg        => undef,
		q_aoh       => [],
		h_aoh       => [],
		dpl_arr     => [],
		max_structs => $opts{max_structs},
	}, $class;
	
	my $num_structs = $self->_ct_read($fn);
	return undef if $num_structs == 0;
	
	$self->_find_all_duplexes();
	
	return $self;
}

###
# PRIVATE METHODS

###
# Test input:
#	my $q_aoh = [
#		{i => 0, let => 'A', state => 'U'               },
#		{i => 1, let => 'G', state => 'D', duplex_i => 7},
#		{i => 2, let => 'G', state => 'U'               },
#		{i => 3, let => 'G', state => 'D', duplex_i => 6},
#		{i => 4, let => 'A', state => 'S'               },
#		{i => 5, let => 'G', state => 'D', duplex_i => 5},
#		{i => 6, let => 'G', state => 'D', duplex_i => 4},
#		{i => 7, let => 'C', state => 'D', duplex_i => 2},
#		{i => 8, let => 'C', state => 'D', duplex_i => 1},
#	];
#
#	my $h_aoh = [
#		{ i => 0, let => 'A', state => 'U'                },
#		{ i => 1, let => 'G', state => 'D', duplex_i => 8 },
#		{ i => 2, let => 'G', state => 'D', duplex_i => 7 },
#		{ i => 3, let => 'A', state => 'S'                },
#		{ i => 4, let => 'C', state => 'D', duplex_i => 6 },
#		{ i => 5, let => 'C', state => 'D', duplex_i => 5 },
#		{ i => 6, let => 'C', state => 'D', duplex_i => 3 },
#		{ i => 7, let => 'C', state => 'D', duplex_i => 1 },
#	];
#
sub _find_all_duplexes
{
	my $self = shift;
	
	my($q_aoh, $h_aoh) = ($self->{q_aoh}, $self->{h_aoh});
	my $dpl_regions = get_dpl_regions($q_aoh, $h_aoh);
	
	foreach my $reg ( @$dpl_regions )
	{
		my @q_reg_aoh   = @$q_aoh[ $reg->{q_start_i}..$reg->{q_end_i} ];
		my @h_reg_aoh_r = reverse @$h_aoh[ $reg->{h_start_i}..$reg->{h_end_i} ];
		
		my($score_mtx, $path_mtx) = do_smith_waterman(\@q_reg_aoh, \@h_reg_aoh_r);
		my($row, $col) = get_start_point(@$score_mtx);
		
		push @{$self->{dpl_arr}}, traceback(\@q_reg_aoh, \@h_reg_aoh_r, $row, $col, $path_mtx);
	}
}

sub print_mtx
{
	my($mtx) = @_;
	
	foreach my $row ( @$mtx )
	{
		print join("\t", @$row)."\n";
	}
}

sub get_dpl_regions
{
	my($q_aoh, $h_aoh) = @_;
	
	my $q_regs   = _aoh_to_duplex_regions( $q_aoh );
	my $h_regs   = _aoh_to_duplex_regions( $h_aoh );
	my $dpl_regs = {};
	foreach my $q_reg ( @$q_regs )
	{
		foreach my $q_i ( $q_reg->{start_i}..$q_reg->{end_i} )
		{
			next if $q_aoh->[$q_i]{state} ne 'D';
			my $h_i = $q_aoh->[$q_i]{duplex_i};
			
			# Look for the h-region that contains the $h_i
			my $h_reg = _get_region_for_position($h_i, $h_regs);
			if( !$h_reg )
			{
				warn "Can't find hit region for query nt: ".Dumper($q_aoh->[$q_i], $h_regs);
				next;
			}
			
			my $key = "$q_reg->{start_i}-$q_reg->{end_i}:$h_reg->{start_i}-$h_reg->{end_i}";
			$dpl_regs->{ $key } = {
				q_start_i => $q_reg->{start_i},
				q_end_i   => $q_reg->{end_i},
				h_start_i => $h_reg->{start_i},
				h_end_i   => $h_reg->{end_i},
			};
		}
	}
	
	return [ sort {$a->{q_start_i} <=> $b->{q_start_i} or $b->{h_start_i} <=> $a->{h_start_i}} values %$dpl_regs ];
}

sub _aoh_to_duplex_regions
{
	my($aoh) = @_;
	
	my $regs = [];
	my $str  = join '', map { $_->{state} } @$aoh;
	$str =~ tr/US/us/;   # for visual purpose
	while( $str =~ /D[DU]*/gi )
	{
		my $prefix = $`;
		my $match  = $&;
		$match =~ s/U+$//i;
		
		my $start_i = length($prefix);
		my $end_i   = $start_i+length($match)-1;
		push @$regs, {start_i => $start_i, end_i => $end_i};
	}
	
	return $regs;
#	return [ sort { $a->{start_i} <=> $b->{start_i} } @$regs ];
}

sub _get_region_for_position
{
	my($pos, $all_regs) = @_;
	
	foreach my $reg ( @$all_regs )
	{
		if( $reg->{start_i} <= $pos && $pos <= $reg->{end_i} )
		{
			return $reg;
		}
	}
	
	return undef;
}

sub traceback
{
	my($q_aoh, $h_aoh_r, $row, $col, $path_mtx) = @_;

	my $ali = {
		q_start_i => undef,
		q_end_i   => $q_aoh->[$row-1]{seq_i},
		h_start_i => $h_aoh_r->[$col-1]{seq_i},
		h_end_i   => undef,
		q_ali     => '',
		h_ali_r   => '',
		match_str => '',
	};
	
	while( $path_mtx->[$row][$col] ne $S )
	{
		die "Something wrong: row = '$row' and col = '$col'" if $row < 0 || $col < 0;
		
		my($q_let, $h_let, $step) = ($q_aoh->[$row-1]{let}, $h_aoh_r->[$col-1]{let}, $path_mtx->[$row][$col]);
		$ali->{q_ali}    .= ($step eq $D || $step eq $U) ? $q_let : '-';
		$ali->{h_ali_r}  .= ($step eq $D || $step eq $L) ? $h_let : '-';
		$ali->{match_str}.=  $step eq $D ? get_match_char($q_let,$h_let) : ' ';
		$ali->{q_start_i} = $q_aoh->[$row-1]{seq_i};
		$ali->{h_end_i}   = $h_aoh_r->[$col-1]{seq_i};
		
		$step eq $D ? ($row-- and $col--) : ($step eq $L ? $col-- : $row--);
	}
	
	$ali->{q_ali}     = reverse $ali->{q_ali};
	$ali->{h_ali_r}   = reverse $ali->{h_ali_r};
	$ali->{match_str} = reverse $ali->{match_str};
	$ali->{len}       = length($ali->{match_str});
	$ali->{n_matches} = $ali->{match_str} =~ s/\|/\|/g;
	$ali->{identity}  = sprintf('%.1f', 100*$ali->{n_matches}/$ali->{len});
	
	return $ali;
}

sub get_match_char
{
	my($q_let, $h_let) = @_;
	
	if( $q_let eq revcomp($h_let) )
	{
		return '|';   # Match
	}
	elsif( ($q_let eq 'G' && $h_let eq 'T') || ($q_let eq 'T' && $h_let eq 'G') )
	{
		return ':';   # Wobble base pair
	}
	else
	{
		 return ' ';
	}
}

sub do_smith_waterman
{
	my($q_aoh, $h_aoh_r) = @_;
	
	my(@score_mtx, @path_mtx);
	$path_mtx[$_]    = [] and $score_mtx[$_]    = [] foreach (0 .. scalar(@$q_aoh));
	$path_mtx[$_][0] = $S and $score_mtx[$_][0] = 0  foreach (0 .. scalar(@$q_aoh));
	$path_mtx[0][$_] = $S and $score_mtx[0][$_] = 0  foreach (0 .. scalar(@$h_aoh_r));
	($path_mtx[0][0], $score_mtx[0][0]) = ($S, 0);
	
	foreach my $row ( 1 .. scalar(@$q_aoh) )
	{
		my $q_nt = $q_aoh->[$row-1];
		foreach my $col ( 1 .. scalar(@$h_aoh_r) )
		{
			my $h_nt = $h_aoh_r->[$col-1];
			_add_value_to_mtx($q_nt, $h_nt, $row, $col, \@score_mtx, \@path_mtx);
		}
	}
	
	return(\@score_mtx, \@path_mtx);
}

sub _add_value_to_mtx
{
	my($q_nt, $h_nt, $row, $col, $score_mtx, $path_mtx) = @_;
	
	if( $q_nt->{state} eq 'D' && $h_nt->{state} eq 'D' )
	{
		# in the $path_mtx we can insert \ or *
		if( $q_nt->{duplex_i} == $h_nt->{seq_i} )
		{
			warn "Something is wrong: ".Dumper($q_nt, $h_nt) if $h_nt->{duplex_i} != $q_nt->{seq_i};
			$score_mtx->[$row][$col] = $score_mtx->[$row-1][$col-1] + 1;
			$path_mtx->[$row][$col]  = $D;
		}
		else
		{
			$score_mtx->[$row][$col] = 0;
			$path_mtx->[$row][$col]  = $S;
		}
	}
	else
	{
		# now we have either D-U, U-D or U-U and in the $path_mtx we can insert -, | or \
		my $diag = $score_mtx->[$row-1][$col-1];
		my $left = $score_mtx->[ $row ][$col-1];
		my $up   = $score_mtx->[$row-1][ $col ];
		$score_mtx->[$row][$col] = (sort { $b <=> $a } $diag, $left, $up, 0)[0];
		if( $score_mtx->[$row][$col] == 0 )
		{
			$path_mtx->[$row][$col] = $S;
		}
		elsif( $score_mtx->[$row][$col] == $diag )
		{
			$path_mtx->[$row][$col] = $D;
		}
		elsif( $score_mtx->[$row][$col] == $left )
		{
			$path_mtx->[$row][$col] = $L;
		}
		else
		{
			$path_mtx->[$row][$col] = $U;
		}
	}
}

sub get_start_point
{
	my(@arr) = @_;
	
	my($max, $max_row, $max_col) = (0,0,0);
	foreach my $row ( 0 .. scalar(@arr)-1 )
	{
		foreach my $col ( 0 .. scalar(@{$arr[0]})-1 )
		{
			if( $arr[$row][$col] > $max )
			{
				$max     = $arr[$row][$col];
				$max_row = $row;
				$max_col = $col;
			}
		}
	}
	
	return($max_row, $max_col);
}

###
# Returns number of structures that have been successfully read
#
sub _ct_read
{
	my $self = shift;
	my($fn) = @_;
	
	open(F, '<', $fn) or die "Can't open file $fn: $!";
	my $all_structs = [];
	while( my $struct = _ct_read_next_struct(\*F) )
	{
		push @$all_structs, $struct;
		last if $self->{max_structs} && scalar(@$all_structs) >= $self->{max_structs};
	}
	close F;
	
	return 0 if scalar(@$all_structs) == 0;
	
	# First, take the data from the Minimal Free Energy (MFE) Structure
	my($q_mfe, $h_mfe) = ($all_structs->[0]{q_aoh}, $all_structs->[0]{h_aoh});
	$self->{enrg}      = $all_structs->[0]{enrg};
	$self->{q_aoh}     = $all_structs->[0]{q_aoh};
	$self->{h_aoh}     = $all_structs->[0]{h_aoh};
	$self->{n_structs} = scalar( @$all_structs );
	
	$self->_find_duplex_nt($all_structs);
	
	$self->_find_unpaired_nt_in('q_aoh', $all_structs);
	$self->_find_unpaired_nt_in('h_aoh', $all_structs);
	
	# verify & set the state
	foreach my $nt ( @{$self->{q_aoh}}, @{$self->{h_aoh}} )
	{
		die "Something is wrong for nt: ".Dumper($nt) if $nt->{duplex} && $nt->{unpaired};
		$nt->{state} = $nt->{duplex} ? 'D' : $nt->{unpaired} ? 'U' : 'S'
	}
	
	return scalar(@$all_structs);
}

sub _find_unpaired_nt_in
{
	my $self = shift;
	my($name, $all_structs) = @_;
	
	# Compute the unpaired probabilities and set the unpaired flags on the corresponding nt
	$_->{unpaired_p} = 0 foreach @{$self->{ $name }};
	$_->{unpaired}   = 0 foreach @{$self->{ $name }};
	for(my $i = 0; $i < scalar( @{$self->{$name}} ); $i++ )
	{
		my $num_structs = 0;   # number of structures where the nt is unpaired
		foreach my $struct ( @$all_structs )
		{
			$num_structs++ if !$struct->{$name}[$i]{paired};
		}
		my $unp_prob = sprintf('%.4f', $num_structs/$self->{n_structs}) + 0;
		$self->{$name}[$i]{unpaired_p} = $unp_prob;
		$self->{$name}[$i]{unpaired}   = 1 if $unp_prob >= 0.9;  # set unpaired flag if we meet the threshold
	}
}

sub _find_duplex_nt
{
	my $self = shift;
	my($all_structs) = @_;
	
	# Compute the duplex probabilities and set the duplex flags on the corresponding nt
	$_->{duplex_p}  = 0 foreach @{$self->{q_aoh}}, @{$self->{h_aoh}};
	$_->{duplex}    = 0 foreach @{$self->{q_aoh}}, @{$self->{h_aoh}};
	for(my $q_i = 0; $q_i < scalar(@{$self->{q_aoh}}); $q_i++ )
	{
		if( my $h_i = $self->{q_aoh}[$q_i]{duplex_i} )
		{
			my $dpl_structs = 0;   # number of structures where the duplex $q_i :: $h_i was present
			foreach my $struct ( @$all_structs )
			{
				$dpl_structs++ if $struct->{q_aoh}[$q_i]{duplex_i} && $struct->{q_aoh}[$q_i]{duplex_i} == $h_i;
			}
			my $dpl_prob = sprintf('%.4f', $dpl_structs/$self->{n_structs}) + 0;
			$self->{h_aoh}[$h_i]{duplex_i} = $q_i;
			$self->{q_aoh}[$q_i]{duplex_p} = $dpl_prob;
			$self->{h_aoh}[$h_i]{duplex_p} = $dpl_prob;
			$self->{q_aoh}[$q_i]{duplex}   = 1 if $dpl_prob == 1;  # set duplex flag only if pairing was observed in 100% of structures
			$self->{h_aoh}[$h_i]{duplex}   = 1 if $dpl_prob == 1;  # set duplex flag only if pairing was observed in 100% of structures
		}
	}
}

sub _ct_read_next_struct
{
	my( $fh ) = @_;
	my $struct = {};
	
	my $head = <$fh>;
	return undef if !$head;
	
	#  4419  ENERGY = -1654.4  49969604_262526738
	my $total_len = undef;    # including the 'I' characters
	if(	$head =~ /(\d+)\s+ENERGY = (\S+)/ )
	{
		($total_len, $struct->{enrg}) = ($1, $2);
	}
	else
	{
		warn("Can't extract energy from the .ct file header line: '$head'");
		return undef;
	}
	
	my $data = read_table_wo_header_from_fh($fh, delim_re => qr/\s+/, max_lines => $total_len);
	
	my $i = 0;
	for(; $i < @$data; $i++)
	{
		my $d = $data->[$i];
		last if $d->[1] eq 'I';
#		push @{$struct->{q_aoh}}, {i => $d->[0], seq_i => $d->[-1], let => $d->[1], paired => $d->[4]};
		push @{$struct->{q_aoh}}, {i => $d->[0]-1, seq_i => $d->[-1]-1, let => $d->[1], paired => $d->[4]};
	}
	
	for(; $i < @$data; $i++)
	{
		my $d = $data->[$i];
		next if $d->[1] eq 'I';
#		push @{$struct->{h_aoh}}, {i => $d->[0], seq_i => $d->[-1], let => $d->[1], paired => $d->[4]};
		push @{$struct->{h_aoh}}, {i => $d->[0]-1, seq_i => $d->[-1]-1, let => $d->[1], paired => $d->[4]};
	}
	
	# Determine inter-molecular duplexes
	my $q_len = scalar(@{$struct->{q_aoh}});
	foreach my $aoh ( @{$struct->{q_aoh}} )
	{
		if( $aoh->{paired} > $q_len )
		{
			$aoh->{duplex_i} = $aoh->{paired} - $q_len - 4;
		}
	}
	
	return $struct;
}

###
# PUBLIC METHODS

sub q_aoh           { return $_[0]->{q_aoh};     }
sub h_aoh           { return $_[0]->{h_aoh};     }
sub get_num_structs { return $_[0]->{n_structs}; }
sub get_full_energy { return $_[0]->{enrg};      }

sub dpl_arr
{
	my $self = shift;
	my(%opts) = @_;
	
	my $q_offset = $opts{q_offset} || 0;
	my $h_offset = $opts{h_offset} || 0;
	warn "Something wrong with the offset: ".Dumper(\%opts) if $q_offset < 0 || $h_offset < 0;
	
	my $dpl_arr  = [];
	foreach my $dpl ( @{$self->{dpl_arr}} )
	{
		my %h = %$dpl;
		
		$h{q_start} = $q_offset + $h{q_start_i} + 1;
		$h{q_end}   = $q_offset + $h{q_end_i}   + 1;
		$h{h_start} = $h_offset + $h{h_start_i} + 1;
		$h{h_end}   = $h_offset + $h{h_end_i}   + 1;
		
		if( $opts{uc_keys} )
		{
			%h = map { uc($_) => $h{$_} } keys %h;
		}
		
		push @$dpl_arr, \%h;
	}
	
	return $dpl_arr;
}
###
#   $arr
#     |
#     |--->{start}  --  1-based
#     |--->{end}
#     |--->{len}
#     |--->{seq}
#
sub q_dpl
{
	my $self = shift;
	
	my $arr = [];
	foreach my $dpl ( @{$self->{dpl_arr}} )
	{
		(my $seq = $dpl->{q_ali}) =~ s/-//g;
		push @$arr, {
			start => $dpl->{q_start_i}+1,
			end   => $dpl->{q_end_i}+1,
			len   => length($seq),
			seq   => $seq,
		};
	}
	
	return $arr;
}

sub h_dpl
{
	my $self = shift;
	
	my $arr = [];
	foreach my $dpl ( @{$self->{dpl_arr}} )
	{
		(my $seq = reverse($dpl->{h_ali_r})) =~ s/-//g;
		push @$arr, {
			start => $dpl->{h_start_i}+1,
			end   => $dpl->{h_end_i}+1,
			len   => length($seq),
			seq   => $seq,
		};
	}
	
	return $arr;
}

sub get_seq
{
	my $self = shift;
	my($name)= shift;
	
	$name   = 'query' if !$name;
	my $key = $name eq 'hit' ? 'h_aoh' : 'q_aoh';
	
	return join('', map { $_->{let} } @{$self->{$key}});
}

sub get_total_dpl_ali_len
{
	my $self = shift;
	
	my $total_len = 0;
	foreach my $dpl ( @{$self->{dpl_arr}} )
	{
		$total_len += length( $dpl->{q_ali} );
	}
	
	return $total_len;
}

sub dpl_len_bp
{
	my $self = shift;

warn "WARNING: the function dpl_len_bp() is depricatred -- consider using get_total_dpl_ali_len()!!!!";
	
	my $len_nt = 0;
	foreach my $h ( @{$self->{q_aoh}}, @{$self->{h_aoh}} )
	{
		$len_nt++ if $h->{duplex};
	}
	warn "WARNING: duplex len in nt is not divisible by 2: $len_nt" if $len_nt % 2 != 0;
	
	return $len_nt/2;
}

sub get_dpl_txt
{
	my $self = shift;
	my(%opts) = @_;
	
	my $q_gene_id = $opts{q_gene_id} || 'Q_GENE';
	my $q_trx_id  = $opts{q_trx_id}  || 'Q_TRX';
	my $q_offset  = $opts{q_offset}  || 0;
	my $h_gene_id = $opts{h_gene_id} || 'H_GENE';
	my $h_trx_id  = $opts{h_trx_id}  || 'H_TRX';
	my $h_offset  = $opts{h_offset}  || 0;
	
	my($txt, $dpl_id) = ('', 'A');
	foreach my $dpl ( @{$self->{dpl_arr}} )
	{
		my $q_start = $q_offset + $dpl->{q_start_i} + 1;
		my $q_end   = $q_offset + $dpl->{q_end_i}   + 1;
		my $h_start = $h_offset + $dpl->{h_start_i} + 1;
		my $h_end   = $h_offset + $dpl->{h_end_i}   + 1;
		$txt .= "== Duplex $dpl_id ==\n";
		$txt .= "Length = $dpl->{len} bp (complementarity=$dpl->{identity}%)\n";
		$txt .= "Query: $q_gene_id [$q_trx_id:$q_start-$q_end]\n";
		$txt .= "5' - $dpl->{q_ali} - 3'\n";
		$txt .= "     $dpl->{match_str}\n";
		$txt .= "3' - $dpl->{h_ali_r} - 5'\n";
		$txt .= "Hit: $h_gene_id [$h_trx_id:$h_start-$h_end]\n";
		$txt .= "\n\n";
		$dpl_id++;
	}
	
	return $txt;
}


###
# SUBROUTINES

1;

