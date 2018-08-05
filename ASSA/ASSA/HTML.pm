package ASSA::HTML;

use strict;
use warnings;

# $Id: HTML.pm 2565 2016-06-13 16:53:23Z antonov $

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
#
# $self
#   |
#   |
#   |--->{verbose}
#   |
#   |--->{install_dir}
#   |--->{html_dir}
#   |--->{css_dir}
#   |--->{img_dir}
#   |
#   |--->{assa}
#   |--->{gpair_aoh}
# 


use Data::Dumper;
use File::Copy qw(copy);

use HTML::Template;
use Bio::Graphics;
use Bio::SeqFeature::Generic;

use ASSA::Lib qw(write_to_file max find_longest_orfs);

###
# CONSTANTS

###
# CONSTRUCTOR
sub new
{
	my $class = shift;
	my($assa, %o) = @_;
	
	my($install_dir)= $INC{"ASSA/ASSA.pm"} =~ /^(.+)\/ASSA\.pm$/;
	
	my $self = bless {
		assa        => $assa,
		gpair_aoh   => $assa->get_gene_pair_aoh,
		install_dir => $install_dir,
		verbose     => $o{verbose},
	}, $class;
}


###
# PRIVATE METHODS

###
# PUBLIC METHODS

sub create_assa_html_dir
{
	my $self = shift;
	my($dir) = @_;
	
	print STDERR "Generating HTML output...\n" if $self->{verbose};
	
	$self->{html_dir} = "$dir/html";
	$self->{img_dir}  = "$dir/images";
	$self->{css_dir}  = "$dir/css";
	foreach my $name ($dir, $self->{html_dir}, $self->{img_dir}, $self->{css_dir})
	{
		mkdir($name) if !-d $name;
	}
	
	copy("$self->{install_dir}/css/style.css", $self->{css_dir});
	
	$self->{IMG_LEGEND_1} = $self->create_legend_1_img();
	$self->{IMG_LEGEND_2} = $self->create_legend_2_img();
	
	my $num = 1;
	foreach my $gpair ( @{$self->{gpair_aoh}} )
	{
		$gpair->{Q_IMG_FN} = "interaction_$num.query.png";
		$gpair->{H_IMG_FN} = "interaction_$num.hit.png";
		$gpair->{HTML_FN}  = "interaction_$num.html";
		
		$self->create_trx_img($gpair, 'Q', $gpair->{Q_IMG_FN}, img_w => 900);
		$self->create_trx_img($gpair, 'H', $gpair->{H_IMG_FN}, img_w => 900);
		$self->create_interaction_page($gpair, $gpair->{HTML_FN});
		$num++;
	}
	
	$self->create_gene_pair_table("index.html");
	
	# Create redirect to the gene pair table
	write_to_file('<meta http-equiv="Refresh" content="0; url=html/index.html">', "$dir/index.html");
}

sub create_interaction_page
{
	my $self = shift;
	my($gpair, $out_fn) = @_;
	
	my $tp = $self->get_template('interaction.html');
	$tp->param(
		q_gene_id    => $gpair->{Q}{GENE_ID},
		q_trx_id     => $gpair->{Q}{TRX_ID},
		h_gene_id    => $gpair->{H}{GENE_ID},
		h_trx_id     => $gpair->{H}{TRX_ID},
		
		img_legend_1 => "../images/$self->{IMG_LEGEND_1}",
		img_legend_2 => "../images/$self->{IMG_LEGEND_2}",
		img_query    => "../images/$gpair->{Q_IMG_FN}",
		img_hit      => "../images/$gpair->{H_IMG_FN}",
		bifold_sites => $gpair->{BIFOLD_SITES},
	);
	$self->create_html_file($tp->output(), "RNA-RNA interaction", $out_fn);
}

sub create_gene_pair_table
{
	my $self = shift;
	my($out_fn) = @_;
	
	my $all_gpairs = [];
	foreach my $gpair ( @{$self->{gpair_aoh}} )
	{
		my($bi_best) = sort { $a->{DDG_PVALUE} <=> $b->{DDG_PVALUE} } @{$gpair->{BIFOLD_SITES}};
		my $num_regs = scalar @{$gpair->{BIFOLD_SITES}};
		push @$all_gpairs, {
			%$gpair,
			Q_GENE_ID  => $gpair->{Q}{GENE_ID},
			H_GENE_ID  => $gpair->{H}{GENE_ID},
			Q_TRX_ID   => $gpair->{Q}{TRX_ID},
			H_TRX_ID   => $gpair->{H}{TRX_ID},
			NUM_REGS   => $num_regs,
			NR_CLASS   => $num_regs > 1 ? 'green' : undef,
			BIFOLD_LEN => $gpair->{DPL_LEN},
			SCORE      => sprintf('%.2f', $gpair->{SCORE}     ),
			DDG_EVALUE => sprintf('%.2e', $gpair->{DDG_EVALUE}),
		};
	}
	
	my $tp = $self->get_template('gene_pair_table.html');
	$tp->param(
		all_gpairs => $all_gpairs,
	);
	$self->create_html_file($tp->output(), "Predicted RNA-RNA interactions", $out_fn);
}

sub get_template
{
	my $self = shift;
	my($tp_name, %opts) = @_;
	
	my $tp = HTML::Template->new(
		filename          => "$self->{install_dir}/html/$tp_name",
		shared_cache      => 0,
		die_on_bad_params => 0,
		loop_context_vars => 1,
		%opts,
	);
	
	return $tp;
}

sub create_html_file
{
	my $self = shift;
	my($text, $title, $out_fn) = @_;
	
	my $tp = $self->get_template('template.html');
	
	$tp->param(
		text   => $text,
		title  => $title,
	);
	
	write_to_file($tp->output, "$self->{html_dir}/$out_fn");
}


sub create_trx_img
{
	my $self = shift;
	my($gpair, $type, $out_fn, %opts) = @_;
	
	my $max_len = max(length($gpair->{Q}{SEQ}), length($gpair->{H}{SEQ}));
	my $trx_len = length($gpair->{$type}{SEQ});
	my $panel = Bio::Graphics::Panel->new(
		-pad_left  => 10,     -length => $max_len,      # to scale the transcript lengths
		-pad_right => 10,     -width  => $opts{img_w},
		-grid      => 1,
	);
	
	if( $type eq 'Q' )
	{
		$self->_trx_img_add_trx_track($panel, $trx_len, direction => 'east');
	}
	elsif( $type eq 'H' )
	{
		$panel->flip(1);
		$self->_trx_img_add_trx_track($panel, $trx_len, direction => 'west');
	}
	
	$self->_trx_img_add_cds_track   ($panel, $gpair->{$type}{SEQ});
	$self->_trx_img_add_alu_track   ($panel, $gpair->{$type}{SEQ_ALU}) if $gpair->{$type}{SEQ_ALU};
	
	my $num = 1;
	foreach my $bi_site ( @{$gpair->{BIFOLD_SITES}} )
	{
		my $blast_site = $self->{assa}->get_blast_site_by_id( $bi_site->{BLAST_SITE_ID} );
		$self->_trx_img_add_reg_track   ($panel, $bi_site, $type, "Region $num");
		$self->_trx_img_add_blast_track ($panel, $blast_site, $type) if $blast_site;
		$self->_trx_img_add_bifold_track($panel, $bi_site, $type);
		$num++;
	}
	
	write_to_file($panel->png, "$self->{img_dir}/$out_fn", binary => 1);
}

sub _trx_img_add_trx_track
{
	my $self = shift;
	my($panel, $trx_len, %o) = @_;
	
	my $feature = Bio::SeqFeature::Generic->new(-start => 1, -end => $trx_len);
	$panel->add_track($feature,
#		-key     => "5' -> $trx->{NAME} ($trx->{EXT_ID}) -> 3'",
		-glyph      => 'arrow',     -tick         => 2,
		-fgcolor    => 'black',     -double       => 0,
		-arrowstyle => 'filled',    $o{direction} => 1,
	);
}

sub _trx_img_add_alu_track
{
	my $self = shift;
	my($panel, $seq_alu) = @_;
	
	my $feature_arr = [];
	while( $seq_alu =~ /N+/g )
	{
		push @$feature_arr, Bio::SeqFeature::Generic->new(
			-start => length($`),
			-end   => length($`) + length($&) - 1,
#			-strand => 1,
		);
	}
	return if @$feature_arr == 0;
	
	$panel->add_track($feature_arr,
		-glyph      => 'generic',
		-stranded   => 0,
		-bgcolor    => 'red',
		-sort_order => 'longest',
	);
}

sub _trx_img_add_cds_track
{
	my $self = shift;
	my($panel, $seq) = @_;
	
	my $cds_arr = find_longest_orfs( $seq );
	return if @$cds_arr == 0 || $cds_arr->[0]{len} < 300;
	
	my $feature = Bio::SeqFeature::Generic->new(
		-start  => $cds_arr->[0]{start},
		-end    => $cds_arr->[0]{len}-1,
		-strand => 1,
	);
	
	$panel->add_track($feature,
		-glyph    => 'generic',
		-bgcolor  => 'green',
		-stranded => 1,
	);
}

sub _trx_img_add_blast_track
{
	my $self = shift;
	my($panel, $site, $type) = @_;
	
	my $feature = Bio::SeqFeature::Generic->new(
		-score        => $site->{IDENTITY},
		-start        => $site->{$type."_START"},
		-end          => $site->{$type."_END"},
	);
	
	my $track = $panel->add_track($feature,
		-glyph     => 'graded_segments',
#		-label     => 1,
		-min_score => 70,
		-max_score => 100,
		-bgcolor   => 'blue',
	);
}

sub _trx_img_add_bifold_track
{
	my $self = shift;
	my($panel, $site, $type) = @_;
	
	my($feature_arr, $dpl_id) = ([], 'A');
	foreach my $dpl ( @{$site->{DPL_ARR}} )
	{
		push @$feature_arr, Bio::SeqFeature::Generic->new(
			-display_name => $dpl_id,
			-start        => $dpl->{$type."_START"},
			-end          => $dpl->{$type."_END"},
		);
		$dpl_id++;
	}
	
	my $track = $panel->add_track($feature_arr,
#		-label     => 1,
		'-glyph'   => 'generic',
		'-bgcolor' => 'orange',
		'-fgcolor' => 'black',
#		'-bump_limit' => 1,      # Maximum number of levels to bump
		'-bump'    => 0,
	);
}

sub _trx_img_add_reg_track
{
	my $self = shift;
	my($panel, $site, $type, $name) = @_;
	
	my $feature = Bio::SeqFeature::Generic->new(
		-display_name => $name,
		-start        => $site->{$type."_START"}+1,
		-end          => $site->{$type."_END"}-1,
	);
	
	$panel->add_track($feature,
		-glyph  => 'anchored_arrow',
		-label  => 1,
	);
}

sub create_legend_1_img
{
	my $self = shift;
	
	my $panel = Bio::Graphics::Panel->new(
		-length => 500, -width  => 400, -pad_left  => 3,
		-key_style => 'between',
	);
	
	my $cds_track = $panel->add_track(
		-key        => 'Coding Sequence (CDS)',
		-glyph      => 'generic',
		-stranded   => 1,
		-bgcolor    => 'green',
	);
	$cds_track->add_feature(Bio::SeqFeature::Generic->new(
		-start  => 1,
		-end    => 100,
		-strand => +1,
	));
	
	my $alu_track = $panel->add_track(
		-key        => 'Alu repeat',
		-glyph      => 'generic',
		-stranded   => 1,
		-bgcolor    => 'red',
	);
	$alu_track->add_feature(Bio::SeqFeature::Generic->new(
		-start  => 1,
		-end    => 100,
		-strand => +1,
	));
	
	my $out_fn = "legend_1.png";
	write_to_file($panel->png, "$self->{img_dir}/$out_fn", binary => 1);
	
	return $out_fn;
}

sub create_legend_2_img
{
	my $self = shift;
	
	my @all_items = (100,90,80,70);
	my $item_len  = 20;
	
	my $panel = Bio::Graphics::Panel->new(
		-length    => scalar(@all_items)*$item_len,
		-width     => 200,
		-pad_left  => 3,
		-key_style => 'between',
	);
	
	my $track = $panel->add_track(
		-key       => 'BLASTn sites (%identity)',
		-label     => 1,
		-glyph     => 'graded_segments',
		-min_score => 70,
		-max_score => 100,     
		-bgcolor   => 'blue',
	);
	
	for(my $i = 0; $i < @all_items; $i++)
	{
		$track->add_feature(Bio::SeqFeature::Generic->new(
			-display_name => $all_items[$i].'%',
			-score        => $all_items[$i],
			-start        => $i * $item_len,
			-end          => $i * $item_len + $item_len - 5,
		));
	}
	
	my $bf_track = $panel->add_track(
		-key        => 'bifold duplex',
		-glyph      => 'generic',
		-bgcolor    => 'orange',
	);
	$bf_track->add_feature(Bio::SeqFeature::Generic->new(
		-start  => 1,
		-end    => 50,
	));
	
	my $out_fn = "legend_2.png";
	write_to_file($panel->png, "$self->{img_dir}/$out_fn", binary => 1);
	
	return $out_fn;
}

1;

