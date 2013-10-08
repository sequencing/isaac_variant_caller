#!/usr/bin/env perl

=head1 LICENSE

Isaac Variant Caller Workflow
Copyright (c) 2009-2013 Illumina, Inc.

This software is provided under the terms and conditions of the
Illumina Open Source Software License 1.

You should have received a copy of the Illumina Open Source
Software License 1 along with this program. If not, see
<https://github.com/sequencing/licenses/>.

=head1 SYNOPSIS

callSmallVariants.pl [options] | --help

=head2 SUMMARY

Run the small variant caller for snvs and indels on a single
chromosome bin.

=cut

use warnings FATAL => 'all';
use strict;

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

if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for strelka_workflow install directory: '$baseDir'");
}
my $libexecDir=File::Spec->catdir($baseDir,'libexec');
my $optDir=File::Spec->catdir($baseDir,'opt');


my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline = join(' ',$0,@ARGV);


my ($chrom, $binId, $configFile);
my $skipHeader = 0;
my $help;

GetOptions( "chrom=s" => \$chrom,
            "bin=s" => \$binId,
            "config=s" => \$configFile,
            "skipHeader" => \$skipHeader,
            "help|h" => \$help) or pod2usage(2);

pod2usage(2) if($help);
pod2usage(2) unless(defined($chrom));
pod2usage(2) unless(defined($binId));
pod2usage(2) unless(defined($configFile));



#
# check all fixed paths (not based on commandline arguments):
#
checkDir($baseDir);
checkDir($libexecDir);
checkDir($optDir);

my $ivcBin=File::Spec->catfile($libexecDir,'starling2');
checkFile($ivcBin,"ivc binary");
my $bgzip9Bin=File::Spec->catfile($optDir,'bgzf_extras','bgzip9');
checkFile($bgzip9Bin,"bgzip9 binary");
my $samtoolsBin = File::Spec->catfile($optDir,'samtools','samtools');
checkFile($samtoolsBin,"samtools binary");



#
# read config and validate values
#
checkFile($configFile,"configuration ini");
my $config  = parseConfigIni($configFile);

for (qw(knownGenomeSize inputBam refFile outDir)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{derived}{$_}));
}
for (qw(depthFilterMultiple minGQX indelMaxRefRepeat minMapq isWriteRealignedBam binSize
        maxInputDepth isSkipDepthFilters extraIvcArguments)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
}

my $outDir = $config->{derived}{outDir};
my $binDir = File::Spec->catdir($outDir,'chromosomes',$chrom,'bins',$binId);
checkDir($outDir,"output");
checkDir($binDir,"output bin");


my $inputBam = $config->{derived}{inputBam};
my $refFile = $config->{derived}{refFile};
checkFile($inputBam,"input BAM");
checkFile($refFile,"reference");


# pull out some config options for convenience:
my $binSize=$config->{user}{binSize};
my $isWriteRealignedBam=$config->{user}{isWriteRealignedBam};
my $knownGenomeSize = $config->{derived}{knownGenomeSize};


my $begin = (int($binId)*$binSize)+1;
my $end = ((int($binId)+1)*$binSize);
#my $end = $begin+100000;  #debug mode


my $useroptions = $config->{user};



#
# setup the command-line:
#
my $base_opts= "--gvcf-file -" .
" --gvcf-max-depth-factor " . $useroptions->{depthFilterMultiple} .
" --gvcf-min-gqx " . $useroptions->{minGQX} .
" --gvcf-max-indel-ref-repeat ". $useroptions->{indelMaxRefRepeat} .
" -bam-file '$inputBam'" .
" -samtools-reference '$refFile'" .
" -bam-seq-name '$chrom'" .
" -report-range-begin $begin -report-range-end $end" .
" -clobber" .
" -min-paired-align-score " . $useroptions->{minMapq} .
" -min-single-align-score " . $useroptions->{minMapq} .
" -bsnp-ssd-no-mismatch 0.35" .
" -bsnp-ssd-one-mismatch 0.6" .
" -min-vexp 0.25" .
" -max-window-mismatch 2 20" .
" -max-indel-size 50" .
" -genome-size $knownGenomeSize";


my $cmd =  "$ivcBin $base_opts";


sub ualignFile($) {
    return File::Spec->catfile($binDir,$_[0] . ".unsorted.realigned.bam");
}
sub alignFile($) {
    return File::Spec->catfile($binDir,$_[0] . ".realigned");
}

my $prefix=getFileBase($inputBam);


if($useroptions->{maxInputDepth} > 0) {
    $cmd .= " --max-input-depth " . $useroptions->{maxInputDepth};
}

if(! $config->{user}{isSkipDepthFilters}) {
    $cmd .= " --chrom-depth-file '" . $config->{derived}{depthFile} . "'";
}

if($isWriteRealignedBam) {
    $cmd .= " -realigned-read-file '" . ualignFile($prefix) . "'";
}

if(defined($useroptions->{extraIvcArguments})){
    my $arg=$useroptions->{extraIvcArguments};
    if($arg !~ /^\s*$/) {
        $cmd .= " " . $arg;
    }
}

if($skipHeader) {
    $cmd .= " --gvcf-skip-header";
}

$cmd .= " |  " . $bgzip9Bin . " -c >| '" . File::Spec->catfile($binDir,"$prefix.genome.vcf.gz") . "'";
$cmd .= " 2>| '" . File::Spec->catfile($binDir,'ivc.stderr') . "'";

executeCmd($cmd,0);


sub sortBam($) {
    my ($label) = @_;
    my $ufile = ualignFile($label);
    if( -f $ufile ) {
        my $afile = alignFile($label);
        my $cmd = "$samtoolsBin sort " . $ufile .  " " . $afile;
        executeCmd($cmd,0);
        unlink($ufile);
    } else {
        logX("Can't find unsorted realigned BAM file: '$ufile'");
    }
}


if($isWriteRealignedBam) {
    sortBam($prefix);
}


1;

__END__

