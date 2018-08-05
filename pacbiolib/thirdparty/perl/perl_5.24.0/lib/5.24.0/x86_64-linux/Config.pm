# This file was created by configpm when Perl was built. Any changes
# made to this file will be lost the next time perl is built.

# for a description of the variables, please have a look at the
# Glossary file, as written in the Porting folder, or use the url:
# http://perl5.git.perl.org/perl.git/blob/HEAD:/Porting/Glossary

package Config;
use strict;
use warnings;
use vars '%Config', '$VERSION';

$VERSION = "5.024000";

# Skip @Config::EXPORT because it only contains %Config, which we special
# case below as it's not a function. @Config::EXPORT won't change in the
# lifetime of Perl 5.
my %Export_Cache = (myconfig => 1, config_sh => 1, config_vars => 1,
		    config_re => 1, compile_date => 1, local_patches => 1,
		    bincompat_options => 1, non_bincompat_options => 1,
		    header_files => 1);

@Config::EXPORT = qw(%Config);
@Config::EXPORT_OK = keys %Export_Cache;

# Need to stub all the functions to make code such as print Config::config_sh
# keep working

sub bincompat_options;
sub compile_date;
sub config_re;
sub config_sh;
sub config_vars;
sub header_files;
sub local_patches;
sub myconfig;
sub non_bincompat_options;

# Define our own import method to avoid pulling in the full Exporter:
sub import {
    shift;
    @_ = @Config::EXPORT unless @_;

    my @funcs = grep $_ ne '%Config', @_;
    my $export_Config = @funcs < @_ ? 1 : 0;

    no strict 'refs';
    my $callpkg = caller(0);
    foreach my $func (@funcs) {
	die qq{"$func" is not exported by the Config module\n}
	    unless $Export_Cache{$func};
	*{$callpkg.'::'.$func} = \&{$func};
    }

    *{"$callpkg\::Config"} = \%Config if $export_Config;
    return;
}

die "$0: Perl lib version (5.24.0) doesn't match executable '$^X' version ($])"
    unless $^V;

$^V eq 5.24.0
    or die sprintf "%s: Perl lib version (5.24.0) doesn't match executable '$^X' version (%vd)", $0, $^V;


sub FETCH {
    my($self, $key) = @_;

    # check for cached value (which may be undef so we use exists not defined)
    return exists $self->{$key} ? $self->{$key} : $self->fetch_string($key);
}

sub TIEHASH {
    bless $_[1], $_[0];
}

sub DESTROY { }

sub AUTOLOAD {
    require 'Config_heavy.pl';
    goto \&launcher unless $Config::AUTOLOAD =~ /launcher$/;
    die "&Config::AUTOLOAD failed on $Config::AUTOLOAD";
}

# tie returns the object, so the value returned to require will be true.
tie %Config, 'Config', {
    archlibexp => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/_output/modulebuilds/thirdparty.src/perl/perl_5.24.0/_output/install/lib/5.24.0/x86_64-linux',
    archname => 'x86_64-linux',
    cc => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/prebuilt.out/3.x/gcc/gcc-4.9.2/libc-2.5/bin/gcc',
    d_readlink => 'define',
    d_symlink => 'define',
    dlext => 'so',
    dlsrc => 'dl_dlopen.xs',
    dont_use_nlink => undef,
    exe_ext => '',
    inc_version_list => ' ',
    intsize => '4',
    ldlibpthname => 'LD_LIBRARY_PATH',
    libpth => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/prebuilt.out/3.x/gcc/gcc-4.9.2/libc-2.5/bin/../lib/gcc/x86_64-libc_2.5-linux-gnu/4.9.2/include-fixed /localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/prebuilt.out/3.x/gcc/gcc-4.9.2/libc-2.5/bin/../lib/gcc/x86_64-libc_2.5-linux-gnu/4.9.2/../../../../x86_64-libc_2.5-linux-gnu/lib /localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/prebuilt.out/3.x/gcc/gcc-4.9.2/libc-2.5/x86_64-libc_2.5-linux-gnu/sysroot/usr/lib /localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/prebuilt.out/3.x/gcc/gcc-4.9.2/libc-2.5/x86_64-libc_2.5-linux-gnu/sysroot/lib',
    osname => 'linux',
    osvers => '3.13.0-24-generic',
    path_sep => ':',
    privlibexp => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/_output/modulebuilds/thirdparty.src/perl/perl_5.24.0/_output/install/lib/5.24.0',
    scriptdir => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/_output/modulebuilds/thirdparty.src/perl/perl_5.24.0/_output/install/bin',
    sitearchexp => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/_output/modulebuilds/thirdparty.src/perl/perl_5.24.0/_output/install/lib/site_perl/5.24.0/x86_64-linux',
    sitelibexp => '/localdisk/build/workspace/S4.0_smrttools-fromsrc_libc-2.5_ubuntu-1404/software/smrtanalysis/_output/modulebuilds/thirdparty.src/perl/perl_5.24.0/_output/install/lib/site_perl/5.24.0',
    so => 'so',
    useithreads => undef,
    usevendorprefix => undef,
    version => '5.24.0',
};
