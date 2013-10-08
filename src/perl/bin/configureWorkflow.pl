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

configureWorkflow.pl --bam FILE --ref FILE --config FILE [options]

This script configures the isaac variant caller workflow for small variant
calling on BAM files. The configuration process will produce an analysis
makefile and directory structure. The makefile can be used to run the
analysis on a workstation or compute cluster via make/qmake or other
makefile compatible process.

=head1 ARGUMENTS

=over 4

=item --bam FILE

Path to indexed input BAM file (required)

=item --ref FILE

Path to indexed reference genome fasta (required)

=item --config FILE

Workflow configuration file. Default config files can be found in
${INSTALL_ROOT}/etc/ for both ELAND and BWA alignments. (required)

=back

=head1 OPTIONS

=over 4

=item --output-dir DIRECTORY

Root of the analysis directory. This script will place all
configuration files in the analysis directory, and after configuration
all results and intermediate files will be written within the analysis
directory during a run. This directory must not already
exist. (default: ./ivcOutput)

=back

=cut

use warnings FATAL => 'all';
use strict;

use Carp;
$SIG{__DIE__} = \&Carp::confess;

use Cwd qw(getcwd);
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


sub usage() { pod2usage(-verbose => 1,
                        -exitval => 2); }

#
# user configuration:
#

my ($inputBam, $refFile, $configFile, $outDir);
my $help;

GetOptions( "bam=s" => \$inputBam,
            "ref=s" => \$refFile,
            "config=s" => \$configFile,
            "output-dir=s" => \$outDir,
            "help|h" => \$help) || usage();

usage() if($help);
usage() unless($argCount);


#
# Validate input conditions:
#

sub checkFileArg($$) {
   my ($file,$label) = @_;

   errorX("Must specify $label file") unless(defined($file));
   checkFile($file,$label);
}

checkFileArg($inputBam,"input BAM");
checkFileArg($refFile,"reference fasta");
checkFileArg($configFile,"configuration ini");

sub makeAbsoluteFilePaths(\$) {
    my ($filePathRef) = @_;

    my ($v,$fileDir,$fileName) = File::Spec->splitpath($$filePathRef);
    if(getAbsPath($fileDir)) {
        errorX("Can't resolve directory path for '$fileDir' from input file argument: '$$filePathRef'");
    }
    $$filePathRef = File::Spec->catfile($fileDir,$fileName);
}

makeAbsoluteFilePaths($inputBam);
makeAbsoluteFilePaths($refFile);
makeAbsoluteFilePaths($configFile);

# also check for BAM index files:
sub checkBamIndex($) {
    my ($file) = @_;
    my $ifile = $file . ".bai";
    if(! -f $ifile) {
        errorX("Can't find index for BAM file '$file'");
    }
}

checkBamIndex($inputBam);


sub checkFaIndex($) {
    my ($file) = @_;
    my $ifile = $file . ".fai";
    if(! -f $ifile) {
        errorX("Can't find index for fasta file '$file'");
    }
    # check that fai file isn't improperly formatted (a la GATK bundle NCBI 37 fai files)
    open(my $FH,"< $ifile") || errorX("Can't open fai file '$ifile'");
    my $lineno=1;
    while(<$FH>) {
          chomp;
          my @F=split();
          if(scalar(@F) != 5) {
              errorX("Unexpected format for line number '$lineno' of fasta index file: '$ifile'\n\tRe-running fasta indexing may fix the issue. To do so, run: \"samtools faidx $file\"");
          }
          $lineno++;
    }
    close($FH);
}

checkFaIndex($refFile);


if(defined($outDir)) {
    if(getAbsPath($outDir)) {
        errorX("Can't resolve path for ouput directory: '$outDir'");
    }
} else {
    $outDir=File::Spec->catdir(Cwd::getcwd(),'ivcOutput');
}

if(-e $outDir) {
    errorX("Output path already exists: '$outDir'");
}

if(getAbsPath($baseDir)) {
    errorX("Can't resolve path for install directory: '$baseDir'");
}

my $samtoolsDir = File::Spec->catdir($optDir,'samtools');
checkDir($libexecDir,"libexec");
checkDir($samtoolsDir,"samtools");

my $callScriptName = "callSmallVariants.pl";
my $finishScriptName = "consolidateResults.pl";
my $callScript = File::Spec->catfile($libexecDir,$callScriptName);
my $finishScript = File::Spec->catfile($libexecDir,$finishScriptName);
my $depthScript = File::Spec->catfile($libexecDir,"getBamAvgChromDepth.pl");
my $countFasta = File::Spec->catfile($libexecDir,"countFastaBases");
my $samtoolsBin = File::Spec->catfile($samtoolsDir,"samtools");
checkFile($callScript,"small variant call script");
checkFile($finishScript,"result consolidation script");
checkFile($depthScript,"bam depth script");
checkFile($countFasta,"fasta scanner");
checkFile($samtoolsBin,"samtools");

#
# Configure bin runs:
#
checkMakeDir($outDir);

#
# Configure bin runs: open and validate config ini
#
my $config = parseConfigIni($configFile);
for my $key (qw(binSize)) {
    unless(exists($config->{user}{$key})) {
        errorX("Config file missing key '$key'");
    }
}

my $binSize = int($config->{user}{binSize});

$config->{derived}{configurationCmdline} = $cmdline;
$config->{derived}{inputBam} = $inputBam;
$config->{derived}{refFile} = $refFile;
$config->{derived}{outDir} = $outDir;

#
# Configure bin runs: check for consistent chrom info between BAMs and reference
#
sub getBamChromInfo($) {
    my $file = shift;
    my $cmd = "$samtoolsBin view -H $file |";
    open(my $FH,$cmd) || errorX("Can't open process $cmd");

    my %info;
    my $n=0;
    while(<$FH>) {
        next unless(/^\@SQ/);
        chomp;
        my @F = split(/\t/);
        scalar(@F) >= 3 || errorX("Unexpected bam header for file '$file'");

        my %h = ();
        foreach (@F) {
            my @vals = split(':');
            $h{$vals[0]} = $vals[1];
        }
        $F[1] = $h{'SN'};
        $F[2] = $h{'LN'};

        my $size = int($F[2]);
        ($size > 0) || errorX("Unexpected chromosome size '$size' in bam header for file '$file'");
        $info{$F[1]}{size} = $size;
        $info{$F[1]}{order} = $n;
        $n++;
    }
    close($FH) || errorX("Can't close process $cmd");
    return %info;
}


my %chromInfo = getBamChromInfo($inputBam);
my @chroms = sort { $chromInfo{$a}{order} <=> $chromInfo{$b}{order} } (keys(%chromInfo));


my %refChromInfo;
{
    logX("Scanning reference genome");
    my $knownGenomeSize=0;
    my $cmd="$countFasta $refFile |";
    open(my $FFH,$cmd) || errorX("Failed to open process '$cmd'");

    while(<$FFH>) {
        chomp;
        my @F = split(/\t/);
        scalar(@F) == 4 || errorX("Unexpected value from '$cmd'");
        $knownGenomeSize += int($F[2]);
        $refChromInfo{$F[1]}{knownSize} = int($F[2]);
        $refChromInfo{$F[1]}{size} = int($F[3]);
    }
    close($FFH) || errorX("Failed to close process '$cmd'");

    #consistency check:
    for my $chrom (@chroms) {
        my $ln = $chromInfo{$chrom}{size};
        my $rln = $refChromInfo{$chrom}{size};
        unless(defined($rln) && ($rln==$ln)) {
            errorX("BAM headers and reference fasta disagree on chromosome: '$chrom'");
        }
        $config->{derived}{"chrom_${chrom}_size"} = $rln;
        $config->{derived}{"chrom_${chrom}_knownSize"} = $refChromInfo{$chrom}{knownSize};
    }
    $config->{derived}{chromOrder} = join("\t",@chroms);

    $config->{derived}{knownGenomeSize} = $knownGenomeSize;
    logX("Scanning reference genome complete");
}


if(! $config->{user}{isSkipDepthFilters}) {

    logX("Estimating chromosome depth");

    my $configDir = File::Spec->catdir($outDir,'config');
    checkMakeDir($configDir);
    $config->{derived}{depthFile} = File::Spec->catdir($configDir,'chrom.depth.txt');
    my $cmd = "perl $depthScript --bam $inputBam >| " . $config->{derived}{depthFile};
    executeCmd($cmd,0);
 
    logX("Estimating chromosome depth complete");
}

#
# Configure bin runs: create directory structure
#
my $resultsDir = File::Spec->catdir($outDir,'results');
checkMakeDir($resultsDir);
if($config->{user}{isWriteRealignedBam}) {
    my $bamDir = File::Spec->catdir($outDir,'realigned');
    checkMakeDir($bamDir);
}
my $chromRootDir = File::Spec->catdir($outDir,'chromosomes');
checkMakeDir($chromRootDir);
for my $chrom (@chroms) {
    my $chromDir = File::Spec->catdir($chromRootDir,$chrom);
    checkMakeDir($chromDir);

    my $chromRef = $chromInfo{$chrom};
    $chromRef->{dir} = $chromDir;
    $chromRef->{binList} = getBinList($chromRef->{size},$binSize);

    my $binRootDir = File::Spec->catdir($chromDir,'bins');
    checkMakeDir($binRootDir);

    for my $binId ( @{$chromRef->{binList}} ) {
        my $binDir = File::Spec->catdir($binRootDir,$binId);
        checkMakeDir($binDir);
    }
}



#
# write run config file:
#
my $runConfigFile;
{
    my $cstr = <<END;
;
; Isaac Variant Caller workflow configuration file
;
; This is an automatically generated file, you probably don't want to edit it. If starting a new run,
; input configuration templates (with comments) can be found in the workflow installation etc/ directory.
;
END

    $cstr .= writeConfigIni($config);

    my $configDir = File::Spec->catdir($outDir,'config');
    checkMakeDir($configDir);
    $runConfigFile = File::Spec->catdir($configDir,'run.config.ini');
    open(my $FH,"> $runConfigFile") || errorX("Can't open file '$runConfigFile'");
    print $FH $cstr;
    close($FH);
}



#
# create makefile
#
my $makeFile = File::Spec->catfile($outDir,"Makefile");
open(my $MAKEFH, "> $makeFile") || errorX("Can't open file: '$makeFile'");

my $completeFile = "task.complete";

print $MAKEFH <<ENDE;
# This makefile was automatically generated by $scriptName
#
# Please do not edit.

script_dir := $libexecDir
call_script := \$(script_dir)/$callScriptName
finish_script := \$(script_dir)/$finishScriptName

config_file := $runConfigFile

analysis_dir := $outDir
results_dir := \$(analysis_dir)/results

ENDE

print $MAKEFH <<'ENDE';

complete_tag := task.complete

finish_task := $(analysis_dir)/$(complete_tag)

get_chrom_dir = $(analysis_dir)/chromosomes/$1
get_bin_task = $(call get_chrom_dir,$1)/bins/$2/$(complete_tag)



all: $(finish_task)
	@$(print_success)


define print_success
echo;\
echo Analysis complete. Final gVCF output can be found in $(results_dir);\
echo
endef


# top level results target:
#
$(finish_task):
	perl $(finish_script) --config=$(config_file) && touch $@


# chromosome targets:
#
ENDE

print $MAKEFH <<ENDE;

# chromosome bin targets:
#
ENDE

my $extraArg = "";

for my $chrom (@chroms) {
    for my $bin (@{$chromInfo{$chrom}{binList}}) {

print $MAKEFH <<ENDE;
chrom_${chrom}_bin_${bin}_task := \$(call get_bin_task,$chrom,$bin)
\$(finish_task): \$(chrom_${chrom}_bin_${bin}_task)
\$(chrom_${chrom}_bin_${bin}_task):
	perl \$(call_script) --config=\$(config_file) --chrom=$chrom --bin=$bin $extraArg&& touch \$@

ENDE

        $extraArg = "--skipHeader ";
    }
}


# If the eval function is available, this is the way we could finish
# the makefile without being so verbose but it doesn't look like qmake
# understands this function.

=cut

print $MAKEFH <<ENDE;

chroms := @chroms

ENDE

for my $chrom (@chroms) {
    print $MAKEFH "${chrom}_bins := " . join(" ",@{$chromInfo{$chrom}{binList}}) . "\n";
}

print $MAKEFH <<'ENDE';

define chrom_task_template
chrom_$1_task := $(call get_chrom_task,$1)
$(finish_task): $$(chrom_$1_task)
$$(chrom_$1_task):
	$$(filter_script) --config=$$(config_file) --chrom=$1 && touch $$@
endef

$(foreach c,$(chroms),$(eval $(call chrom_task_template,$c)))


# chromosome bin targets:
#
define chrom_bin_task_template
chrom_$1_bin_$2_task := $(call get_bin_task,$1,$2)
$$(chrom_$1_task): $$(chrom_$1_bin_$2_task)
$$(chrom_$1_bin_$2_task):
	$$(call_script) --config=$$(config_file) --chrom=$1 --bin=$2 && touch $$@
endef

$(foreach c,$(chroms), \
    $(foreach b,$($c_bins),$(eval $(call chrom_bin_task_template,$c,$b))) \
 )

ENDE

=cut



print <<END;


Successfully configured analysis and created makefile '$makeFile'.

To run the analysis locally using make, run:

make -C $outDir

...or:

cd $outDir
make

END

1;

__END__

