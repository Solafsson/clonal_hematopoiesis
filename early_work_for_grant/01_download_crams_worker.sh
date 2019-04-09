#!/bin/bash

output_dir=$1
jobFile=$2
samtools_path=$3
script_dir=$4
clusterOutput=$5
referenceFile=$6
jobIndex=$7

indexSamplePath=$( awk -v jI="$jobIndex" 'NR==jI {print $2; exit}' < ${jobFile} )
indexSample=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${jobFile} )
indexEGANname=$( awk -v jI="$jobIndex" 'NR==jI {print $4; exit}' < ${jobFile} )

echo "iget -K -n $(ils -l ${indexSamplePath}/${indexSample}.cram | awk '/green/ {print $2; exit}') ${indexSamplePath}/${indexSample}.cram ${output_dir}${indexEGANname}.cram"

retriesfile=$(mktemp so11_irodsXXX)
if iget --retries 3 -X ${retriesfile} -T -f -K  --lfrestart ${retriesfile} -n  $(ils -l ${indexSamplePath}/${indexSample}.cram | awk '/green/ {print $2; exit}') ${indexSamplePath}/${indexSample}.cram ${output_dir}${indexEGANname}.cram; then

    if [[ -s ${output_dir}${indexEGANname}.cram ]]; then
        #${samtools_path} index ${output_dir}${indexEGANname}.cram ${output_dir}${indexEGANname}.cram.crai
        echo "Would have indexed"
    else
        echo "File ${output_dir}${indexEGANname}.cram  is 0 length. Downloading the other copy"
        if iget --retries 3 -X ${retriesfile} -T -f -K  --lfrestart ${retriesfile} -n  $(ils -l ${indexSamplePath}/${indexSample}.cram | awk '/red/ {print $2; exit}') ${indexSamplePath}/${indexSample}.cram ${output_dir}${indexEGANname}.cram; then
            #${samtools_path} index ${output_dir}${indexEGANname}.cram ${output_dir}${indexEGANname}.cram.crai
            echo "Would have indexed"
        else
            echo "Problem with the other copy as well"
            exit 1
        fi
    fi
else
    echo "iget failed for ${indexSamplePath}/${indexSample}.cram"
fi

ls -lah ${output_dir}${indexEGANname}.cram
rm -f ${retriesfile}

echo "Submitting job:"
echo "bsub -J 'cram2bam' -M250 -R'span[hosts=1] select[mem>250] rusage[mem=250]' \
-e ${clusterOutput}cram2bam_errors.%J  -o ${clusterOutput}cram2bam_output.%J  \
bash ${script_dir}02_convert_cram2bam.sh ${output_dir} ${output_dir}${indexEGANname}.cram ${samtools_path} ${referenceFile} ${script_dir}"


bsub -J "cram2bam" -M250 -R'span[hosts=1] select[mem>250] rusage[mem=250]' \
-e ${clusterOutput}cram2bam_errors.%J  -o ${clusterOutput}cram2bam_output.%J -n8 \
bash ${script_dir}02_convert_cram2bam.sh ${output_dir} ${output_dir}${indexEGANname}.cram ${samtools_path} ${referenceFile} ${script_dir}


exit $?
