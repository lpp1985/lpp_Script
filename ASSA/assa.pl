#!/usr/bin/perl --

use strict;
use warnings;

my $SVN_STR  = '$Id: assa.pl 2581 2016-09-13 14:20:59Z antonov $';
my $ASSA_VER = 0.02;

###
# Ivan Antonov (antonov1986@gmail.com)
#

$| = 1; # Turn off buffering

use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw(remove_tree);
use Cwd qw(abs_path);
use IPC::Open2;

use ASSA::ASSA;

###
# CONSTANTS
my $GROUP_BY_GENE   = 0;
my $ALLOW_OVERLAPS  = 1;
my $DISABLE_BBO     = 1;
my $FILTER_ISOFORMS = 0;
my $FILTER_ALU      = 0;
my $FILTER_CIS      = 0;
my $DB_SIZE         = 1;
my $EVALUE_DDG      = 1;
my $PVALUE_DDG      = 1;
my $PVALUE_LEN      = 1;
my $PVALUE_PCC      = 1e-6;
my $HTML_DIR        = undef;
my $FLANK_MIN       = 25;
my $FLANK_MAX       = 50;
my $FLANK_STEP      = 25;
my $SEED_LEN        = 10;
my $BLAST_SITE_LEN  = undef;
my $NUM_STRUCTS     = 1;
my $NUM_THREADS     = 1;
my $E_PROFILES      = undef;
my $Q_ALU_FN        = undef;
my $H_ALU_FN        = undef;
my $OUT_DIR         = undef;
my $TMP_DIR         = '.';
my $SILENT          = 0;

###
# Parse input data
GetOptions(
	'group_by_gene'   => \$GROUP_BY_GENE,
	'allow_overlaps'  => \$ALLOW_OVERLAPS,
	'disable_bbo'     => \$DISABLE_BBO,
	'filter_isoforms' => \$FILTER_ISOFORMS,
	'filter_alu'      => \$FILTER_ALU,
	'filter_cis'      => \$FILTER_CIS,
	'db_size=i'       => \$DB_SIZE,
	'evalue=s'        => \$EVALUE_DDG,
	'pvalue_ddg=s'    => \$PVALUE_DDG,
	'pvalue_len=s'    => \$PVALUE_LEN,
	'pvalue_pcc=s'    => \$PVALUE_PCC,
	'html=s'          => \$HTML_DIR,
	'seed_len=i'      => \$SEED_LEN,
	'blast_len=i'     => \$BLAST_SITE_LEN,
	'flank_min=i'     => \$FLANK_MIN,
	'flank_max=i'     => \$FLANK_MAX,
	'flank_step=i'    => \$FLANK_STEP,
	'num_structs=i'   => \$NUM_STRUCTS,
	'num_threads=i'   => \$NUM_THREADS,
	'q_alu_fn=s'      => \$Q_ALU_FN,
	't_alu_fn=s'      => \$H_ALU_FN,
	'e_profiles=s'    => \$E_PROFILES,
	'out_dir=s'       => \$OUT_DIR,
	'tmp_dir=s'       => \$TMP_DIR,
	'silent'          => \$SILENT,
) || die usage();

die usage() if @ARGV!=2;
die usage("ERROR: directory '$E_PROFILES' doesn't exist!\n") if $E_PROFILES && !-d $E_PROFILES;

###
my $START_TIME = time;

run(
	q_fn            => $ARGV[0],
	h_fn            => $ARGV[1],
	group_by_gene   => $GROUP_BY_GENE,
	filter_bbo      => $DISABLE_BBO    ? 0 : 1,
	filter_overlaps => $ALLOW_OVERLAPS ? 0 : 1,
	filter_isoforms => $FILTER_ISOFORMS,
	filter_alu      => $FILTER_ALU,
	filter_cis      => $FILTER_CIS,
	db_size         => $DB_SIZE,
	evalue_ddg      => $EVALUE_DDG,
	pvalue_ddg      => $PVALUE_DDG,
	pvalue_len      => $PVALUE_LEN,
	pvalue_pcc      => $PVALUE_PCC,
	q_alu_fn        => $Q_ALU_FN,
	h_alu_fn        => $H_ALU_FN,
	flank_min       => $FLANK_MIN,
	flank_max       => $FLANK_MAX,
	flank_step      => $FLANK_STEP,
	seed_len        => $SEED_LEN,
	min_site_len    => $BLAST_SITE_LEN,
	num_structs     => $NUM_STRUCTS,
	num_threads     => $NUM_THREADS,
	e_profiles      => $E_PROFILES,
	html_dir        => $HTML_DIR,
	out_dir         => $OUT_DIR,
	tmp_root        => $TMP_DIR,
	verbose         => $SILENT ? 0 : 1,
);

warn "\nElapsed time: ".(time-$START_TIME)." sec\n" if !$SILENT;
###

###
# SUBROUTINES
sub run
{
	my %opts = @_;
	
	$opts{tmp_dir}  = create_tmp_dir(root_dir => $opts{tmp_root}, prefix => '__assa_');
	$opts{q_alu_fn} = run_repeat_masker($opts{q_fn},%opts) if $opts{filter_alu} && !$opts{q_alu_fn};
	$opts{h_alu_fn} = run_repeat_masker($opts{h_fn},%opts) if $opts{filter_alu} && !$opts{h_alu_fn};
	
	my $assa = ASSA::ASSA->new(%opts);
	
	# Run BLAST & apply FILTER 1 (100% interacting isoforms)
	$assa->run_blast($opts{seed_len}, %opts);
	
	# FILTER 2: gene pairs where both genes have Alu in any direction
	$assa->remove_gene_pairs_with_alu() if $opts{filter_alu};
	
	# FILTER 3: filter out cis-gene pairs
	$assa->remove_cis_gene_pairs() if $opts{filter_cis};
	
	# FILTER 4: keep only gene pairs with correlated expression profiles
	$assa->filter_gene_pairs_by_coexpression($opts{e_profiles}, $opts{pvalue_pcc}) if $opts{e_profiles};
	
	for(my $len = $opts{flank_min}; $len <= $opts{flank_max}; $len += $opts{flank_step})
	{
		$assa->run_bifold($len, "flank_$len",  %opts);
	}
	
	$opts{group_by_gene} ? output_gene_table($assa, %opts) : output_site_table($assa, %opts);
	
	if( $opts{html_dir} )
	{
		require ASSA::HTML;
		ASSA::HTML->new($assa,%opts)->create_assa_html_dir($opts{html_dir});
	}
	
	remove_tree($opts{tmp_dir}) or warn "Couldn't remove tmp dir $opts{tmp_dir}";
}

sub output_gene_table
{
	my($assa, %opts) = @_;
	
	my @head  = qw(Q_NAME    H_NAME    Q_TRX_LEN    H_TRX_LEN                                        );
	push @head, qw(PCC       PCC_NUM   PCC_PVALUE                                                    ) if $assa->has_coexpression;
#	push @head, qw(Q_START   Q_END     H_START      H_END     DPL_LEN   DDG   DDG_PVALUE   DDG_EVALUE) if $assa->has_bifold;
	push @head, qw(NUM_SITES DPL_LEN   SCORE        DDG_EVALUE                                       ) if $opts{group_by_gene};
	
	my $site_aoh = $opts{group_by_gene} ? $assa->get_gene_pair_aoh : $assa->get_bifold_site_aoh;
	
	print join("\t", @head)."\n";
	foreach my $site ( @$site_aoh )
	{
		my $coexpression = $site->{COEXPRESSION};
		
		my %h = (
			Q_NAME         => $site->{Q}{TRX_ID},
			H_NAME         => $site->{H}{TRX_ID},
			Q_TRX_LEN      => length($site->{Q}{SEQ}),
			H_TRX_LEN      => length($site->{H}{SEQ}),
			
			PCC            => $coexpression->{PCC},
			PCC_NUM        => $coexpression->{PCC_NUM},
			PCC_PVALUE     => $coexpression->{PCC_PVALUE},
			
			NUM_SITES      => $site->{NUM_SITES},
			DPL_LEN        => $site->{DPL_LEN},
			SCORE          => $site->{SCORE},
			DDG_EVALUE     => $site->{DDG_EVALUE},
			
#			Q_START        => $site->{Q_START},
#			Q_END          => $site->{Q_END},
#			H_START        => $site->{H_START},
#			H_END          => $site->{H_END},
#			DDG            => $site->{DDG},
#			DDG_PVALUE     => $site->{DDG_PVALUE},
		);
		
		print join("\t", map { defined $h{$_} ? $h{$_} : '-' } @head)."\n";
	}
}

sub output_site_table
{
	my($assa, %opts) = @_;
	
	my @head  = qw(Q_NAME    H_NAME    Q_TRX_LEN    H_TRX_LEN                           );
	push @head, qw(Q_START   Q_END     H_START      H_END     DPL_LEN   DDG   DDG_PVALUE) if $assa->has_bifold;
	
	my $site_aoh = $assa->get_bifold_site_aoh;
	
	print join("\t", @head)."\n";
	foreach my $site ( @$site_aoh )
	{
		my %h = (
			Q_NAME         => $site->{Q}{TRX_ID},
			H_NAME         => $site->{H}{TRX_ID},
			Q_TRX_LEN      => length($site->{Q}{SEQ}),
			H_TRX_LEN      => length($site->{H}{SEQ}),
			
			DPL_LEN        => $site->{DPL_LEN},
			
			Q_START        => $site->{Q_START},
			Q_END          => $site->{Q_END},
			H_START        => $site->{H_START},
			H_END          => $site->{H_END},
			DDG            => $site->{DDG},
			DDG_PVALUE     => $site->{DDG_PVALUE},
		);
		
		print join("\t", map { defined $h{$_} ? $h{$_} : '-' } @head)."\n";
	}
}

sub run_repeat_masker
{
	my($fasta_fn, %opts) = @_;
	
	my $masked_fn = "$opts{tmp_dir}/$fasta_fn.masked";
	my $rm_str = "RepeatMasker -qq -alu -species human -dir $opts{tmp_dir} $fasta_fn";
	$opts{verbose} ? system("$rm_str 1>&2") : `$rm_str`;
	
	$masked_fn = undef if !-e $masked_fn;   # No Alu repeats were found!
	
	return $masked_fn;
}

sub create_tmp_dir
{
	my(%opts) = @_;
	$opts{root_dir} ||= '.';
	$opts{prefix}   ||= '';
	
	my $dir;
	while(1)
	{
		$dir = $opts{root_dir}.'/'.$opts{prefix}.int(rand(1000000));
		last unless -d $dir;
	}
	
	mkdir($dir);
	system("chmod 0777 $dir");
	
	return abs_path($dir);
}

sub usage
{
	my($msg) = @_;
	$msg = $msg ? $msg."\n" : '';
	
	# $Id: assa.pl 2581 2016-09-13 14:20:59Z antonov $
	my($revision, $date) = $SVN_STR =~ /assa.pl\s+(\S+)\s+(\S+)/;
	
	my $script = File::Spec->splitpath($0);

#    --pvalue_len  <NUM>  --  total duplex length p-value threshold. Default value: $PVALUE_LEN
	return"$msg
ASSA, version $ASSA_VER, revision $revision ($date)

USAGE:
    $script   [OPTIONS]   <QUERIES.fna>   <TARGETS.fna>   >   <OUTPUT.txt>

OPTIONS:
    --group_by_gene      --  for each gene pair output the best site only

    --disable_bbo        --  disable the bbo (BLAST-bifold overlap) filter
    --allow_overlaps     --  do not filter overlapping BLASTn sites/bifold regions
    --filter_isoforms    --  100% of all isoforms of two genes must have the interactions.
                             Fasta info-line must begin with: '>gene_id|XXX|...'
    --filter_cis         --  filter out cis-interactions
    --filter_alu         --  run RepeatMasker & filter out possible Alu-based interactions
    --q_alu_fn <FN.fna>  --  transcript sequences from <QUERIES.fna> with already masked Alu repeats
    --t_alu_fn <FN.fna>  --  transcript sequences from <TARGETS.fna> with already masked Alu repeats

    --num_threads <INT>  --  number of threads to use. Default value: $NUM_THREADS
    --num_structs <INT>  --  number of bifold MFE structures to consider. Default value: $NUM_STRUCTS
    --db_size     <INT>  --  effective number of sequences in the target database.
    --seed_len    <INT>  --  seed length (-word_size) for the BLASTn search. Default value: $SEED_LEN
    --blast_len   <INT>  --  min length of BLASTn site. By default it is computed based on site's identity and GC%.
    --flank_min   <INT>  --  length of the shortest flank for gradual filtering by bifold. Default value: $FLANK_MIN
    --flank_max   <INT>  --  length of the longest flank for gradual filtering by bifold. Default value: $FLANK_MAX
    --flank_step  <INT>  --  step for the gradual increasement of flank length in bifold filtering. Default value: $FLANK_STEP
    --evalue      <NUM>  --  E-value threshold. Default value: $EVALUE_DDG
    --pvalue_ddg  <NUM>  --  interaction energy p-value threshold. Default value: $PVALUE_DDG
    --pvalue_pcc  <NUM>  --  Pearson correlation p-value threshold. Default value: $PVALUE_PCC
    --e_profiles  <DIR>  --  folder with expression profiles, file names should be equal to gene names
                             from the <QUERIES.fna> and <TARGETS.fna>; each file contains a single column
                             with numbers (no header, the same number of lines in each file)
    
    --html        <DIR>  --  generate the HTML-output
    --out_dir     <DIR>  --  save files with detailed info to this folder
    --tmp_dir     <DIR>  --  where the temporary folder will be created. Default value: $TMP_DIR
    --silent

";
}

