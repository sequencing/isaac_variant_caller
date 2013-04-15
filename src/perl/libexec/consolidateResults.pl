#!/usr/bin/env perl

=head1 LICENSE

Isaac Variant Caller Workflow
Copyright (c) 2009-2013 Illumina, Inc.

This software is provided under the terms and conditions of the
Illumina Open Source Software License 1.

You should have received a copy of the Illumina Open Source
Software License 1 along with this program. If not, see
<https://github.com/downloads/sequencing/licenses/>.

=head1 SYNOPSIS

consolidateVariants.pl [options] | --help

=head2 SUMMARY

Aggregate final results from all chromosomes

=cut

use warnings FATAL => 'all';
use strict;

use Carp;
$SIG{__DIE__} = \&Carp::confess;

use File::Spec;
use File::Temp;
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
my $optDir=File::Spec->catdir($baseDir,'opt');


my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline=join(' ',$0,@ARGV);


my $configFile;
my $help;

GetOptions( "config=s" => \$configFile,
            "help|h" => \$help) or pod2usage(2);

pod2usage(2) if($help);
pod2usage(2) unless(defined($configFile));

#
# check fixed paths
#
my $samtoolsBin = File::Spec->catfile($optDir,'samtools','samtools');
my $tabixBin = File::Spec->catfile($optDir,'tabix','tabix');
my $bcatBin = File::Spec->catfile($optDir,'bgzf_extras','bgzf_cat');
checkFile($samtoolsBin,"samtools binary");
checkFile($tabixBin,"tabix binary");
checkFile($bcatBin,"bgzf_cat binary");


#
# read config and validate values
#
checkFile($configFile,"configuration ini");
my $config  = parseConfigIni($configFile);


for (qw(inputBam outDir chromOrder)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{derived}{$_}));
}
for (qw(isWriteRealignedBam binSize)) {
    errorX("Undefined configuration option: '$_'") unless(defined($config->{user}{$_}));
}

my $userconfig = $config->{user};

my @chromOrder = split(/\t/,$config->{derived}{chromOrder});
for my $chrom (@chromOrder) {
    my $chromSizeKey = "chrom_" . $chrom . "_size";
    errorX("Undefined configuration option: '$_'") unless(defined($chromSizeKey));
}

# check that all input directories exist:
my $outDir = $config->{derived}{outDir};
checkDir($outDir,"output");
for my $chrom (@chromOrder) {
    my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);
    checkDir($chromDir,"input chromosome");

    my $chromSizeKey = "chrom_" . $chrom . "_size";
    my $binList = getBinList($config->{derived}{$chromSizeKey},$userconfig->{binSize});
    for my $binId (@$binList) {
        my $dir = File::Spec->catdir($chromDir,'bins',$binId);
        checkDir($dir,"input bin");
    }
}



#
# Returns a reference to a list of bin level filenames, in expected chromosome/bin order:
#
sub getBinFileList($$) {
    my ($fileName,$label) = @_;

    my @fileList;
    for my $chrom (@chromOrder) {
        my $chromDir = File::Spec->catdir($outDir,'chromosomes',$chrom);

        my $chromSizeKey = "chrom_" . $chrom . "_size";
        my $binList = getBinList($config->{derived}{$chromSizeKey},$userconfig->{binSize});
        for my $binId (@$binList) {
            my $binDir = File::Spec->catdir($chromDir,'bins',$binId);
            my $path = File::Spec->catfile($binDir,$fileName);
            checkFile($path,"bin-level $label");

            push @fileList,$path;
        }
    }

    return \@fileList;
}



# suffix used for large result file intermediates:
my $itag = ".incomplete";



#
# consolidate all bin-level bam files into a single realigned bam output:
#
sub consolidateBam($) {
    my ($fileName) = @_;

    my $reDir = File::Spec->catdir($outDir,'realigned');
    checkMakeDir($reDir);

    my $bamListRef = getBinFileList($fileName,"realigned bam");
    return unless(scalar(@{$bamListRef}));

    my $headerFH = File::Temp->new();
    my $getHeaderCmd = "bash -c '$samtoolsBin view -H ".$bamListRef->[0]." > $headerFH'";
    executeCmd($getHeaderCmd);

    my $allFile = File::Spec->catfile($reDir,$fileName . $itag);
    my $cmd="$samtoolsBin merge -h $headerFH $allFile ". join(" ",@{$bamListRef});
    executeCmd($cmd);

    my $allFileFinished = File::Spec->catfile($reDir,$fileName);
    checkMove($allFile,$allFileFinished);

    my $indexCmd="$samtoolsBin index $allFileFinished";
    executeCmd($indexCmd);

    unlink(@{$bamListRef});
}



#
# consolidate all bin-level vcf files into a single sorted vcf output:
#
sub consolidateVcf($) {
    my ($fileName) = @_;

    my $vcfListRef = getBinFileList($fileName,"genome vcf");
    unless(scalar(@{$vcfListRef}) > 0) {
        errorX("No bin-level gVCF output");
    }

    my $vcfFinished = File::Spec->catfile($outDir,'results',$fileName);
    if(scalar(@{$vcfListRef}) > 1) {
        my $bcatCmd = "bash -c '$bcatBin -o \'$vcfFinished\'";
        for my $file (@{$vcfListRef}) { $bcatCmd .= " \'$file\'"; }
        $bcatCmd .= "'";
        executeCmd($bcatCmd);
    } else {
        checkMove($vcfListRef->[0],$vcfFinished);
    }

    my $indexCmd="$tabixBin -p vcf '$vcfFinished'";
    executeCmd($indexCmd);

    unlink(@{$vcfListRef});
}



my $prefix = getFileBase($config->{derived}{inputBam});
my $outputVcfName = "$prefix.genome.vcf.gz";
consolidateVcf($outputVcfName);


if($userconfig->{isWriteRealignedBam}) {
    consolidateBam("$prefix.realigned.bam");
}



1;

__END__

