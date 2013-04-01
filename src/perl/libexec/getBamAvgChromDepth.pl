#!/usr/bin/env perl

=head1 LICENSE

Isaac Variant Caller Software
Copyright (c) 2009-2012 Illumina, Inc.

This software is provided under the terms and conditions of the
Illumina Open Source Software License 1.

You should have received a copy of the Illumina Open Source
Software License 1 along with this program. If not, see
<https://github.com/downloads/sequencing/licenses/>.

The distribution includes the code libraries listed below in the
'redist' sub-directory. These are distributed according to the
licensing terms governing each library.

=head1 SYNOPSIS

getBamAvgChromDepth.pl [options] | --help

Estimate average chromosome depth from a WGS BAM file. This will
not work correctly for exome or other targeted data.

=head1 ARGUMENTS

=over 4

=item --bam FILE

Path to indexed input BAM file (required)

=cut


use warnings FATAL => 'all';
use strict;

# debug with stacktrace:
use Carp;
$SIG{__DIE__} = \&Carp::confess;

use File::Spec;
use Getopt::Long;
use Pod::Usage;


my $baseDir;
my $libDir;
BEGIN {
    my $thisDir=(File::Spec->splitpath($0))[1];
    $baseDir=File::Spec->catdir($thisDir,File::Spec->updir());
    $libDir=File::Spec->catdir($baseDir,'lib');
}
use lib $libDir;
use Utils;

#
# check fixed dependencies:
#
if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for strelka_workflow install directory: '$baseDir'");
}
my $optDir=File::Spec->catdir($baseDir,'opt');
checkDir($optDir);

my $samtoolsBin = File::Spec->catfile($optDir,'samtools','samtools');
checkFile($samtoolsBin,"samtools binary");

sub usage() { pod2usage(-verbose => 1,
                        -exitval => 2); }

my $bamfile;
my $help;

GetOptions(
    "bam=s" => \$bamfile,
    "help|h" => \$help) or usage();

usage() if($help);
usage() unless(defined($bamfile));

checkFile($bamfile,"bam");


my $cmd1="$samtoolsBin idxstats $bamfile";
open(my $FP1,"$cmd1 |");

my %chrom;
my @chroms;

while(<$FP1>) {
    my @X =split(/\t/);
    next if($X[0] eq "*");
    $chrom{$X[0]} = [ $X[1] , $X[2] ];
    push @chroms, $X[0];
}

close($FP1);

# pass 0 is a subsampled approximation of read length,
# if that fails, then pass 1 runs the exact computation
#
my $length = 0;
my $count = 0;
for my $pass ((0,1)) {
    my $cmd2;
    if($pass == 0) {
        $cmd2="$samtoolsBin view -F 4 -s 0.1 $bamfile";
    } else {
        $cmd2="$samtoolsBin view -F 4 $bamfile";
    }

    my $sid=open(my $FP2,"$cmd2 |");

    $length = 0;
    $count = 0;
    while(<$FP2>) {
        my @F = split(/\t/,$_,7);
        next unless($F[5] =~ /[0-9]+M/);
        $length += $1 while($F[5] =~ /([0-9]+)M/g);
        $count++;
        last if($count >= 200000);
    }
    kill('INT',$sid);
    close($FP2);

    last if($count > 100000);
    logX("Poor read length approximation results. Count: '$count' Rerunning exact estimate"); # rerun pass 1 for exact read length
}


my $avg_length = ($length/$count);

for my $c (@chroms) {
    next if($chrom{$c}->[0] < $avg_length);
    my $depth = $chrom{$c}->[1]*$avg_length / $chrom{$c}->[0];
    printf "%s\t%.3f\t%s\t%.3f\n",$c,$depth,$count,$avg_length;
}


1;

__END__

