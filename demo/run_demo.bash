#!/usr/bin/env bash

#
# Execute small ivc demonstration/verification run
#

set -o nounset
set -o pipefail

script_dir=$(dirname $0)
root_dir=$script_dir/..
data_dir=$script_dir/data
expected_dir=$script_dir/expected_results

analysis_dir=./ivcDemoAnalysis

config_script=$root_dir/configureWorkflow.pl

demo_config_file=$script_dir/ivc_demo_config.ini





if [ ! -e $config_script ]; then
    cat<<END 1>&2

ERROR: The ivc workflow must be built prior to running demo.
       See installation instructions in '$root_dir/README'

END
    exit 2
fi


#
# Step 1: configure demo
#
if [ -e $analysis_dir ]; then
    cat<<END 1>&2

ERROR: Demo analysis directory already exists: '$analysis_dir' 
       Please remove/move this to rerun demo.

END
    exit 2
fi

cmd="$config_script \
--bam=$data_dir/NA12891_dupmark_chr20_region.bam \
--ref=$data_dir/chr20_860k_only.fa \
--config=$demo_config_file \
--output-dir=$analysis_dir"

echo 1>&2
echo "**** Starting demo configuration." 1>&2
echo "**** Configuration cmd: '$cmd'" 1>&2
echo 1>&2
$cmd

if [ $? -ne 0 ]; then
    echo 1>&2
    echo "ERROR: Demo configuration step failed" 1>&2
    echo 1>&2
    exit 1
else
    echo 1>&2
    echo "**** Completed demo configuration." 1>&2
    echo 1>&2
fi


#
# Step 2: run demo (on single local core):
#
#stderrlog=$analysis_dir/make.stderr.log
cmd="make -C $analysis_dir"
echo 1>&2
echo "**** Starting demo workflow execution." 1>&2
echo "**** Workflow cmd: '$cmd'" 1>&2
echo 1>&2
$cmd


if [ $? -ne 0 ]; then
    cat<<END 1>&2

ERROR: Workflow execution step failed

END
#        See make error log file: '$stderrlog'
    exit 1
else
    echo 1>&2
    echo "**** Completed demo workflow execution." 1>&2
    echo 1>&2
fi


#
# Step 3: Compare results to expected calls
#
results_dir=$analysis_dir/results
echo 1>&2
echo "**** Starting comparison to expected results." 1>&2
echo "**** Expected results dir: $expected_dir" 1>&2
echo "**** Demo results dir: $results_dir" 1>&2
echo 1>&2

filter_variable_metadata() {
    awk '!/^##(fileDate|startTime|cmdline|source_version|reference)/'
}

compressed_vcf_compare() {
    gzip -dc | filter_variable_metadata
}

for f in $(ls $expected_dir/*.gz); do
    file=$(basename $f)
    efile=$expected_dir/$file
    rfile=$results_dir/$file
    diff <(compressed_vcf_compare < $efile) <(compressed_vcf_compare < $rfile)

    if [ $? -ne 0 ]; then
        cat<<END 1>&2

ERROR: Found difference between demo and expected results in file '$f'.
       Expected file: $efile
       Demo results file: $rfile

END
        exit 1
    fi
done

echo 1>&2
echo "**** No differences between expected and computed results." 1>&2
echo 1>&2

echo 1>&2
echo "**** Demo/verification successfully completed" 1>&2
echo 1>&2

