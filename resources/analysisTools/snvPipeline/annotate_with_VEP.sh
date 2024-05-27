#!/bin/bash

## Annotate the variants with VEP

## To run
## LOCAL: sh annotate_variant_with_VEP.sh snvs_{pid}_somatic_snvs_conf_8_to_10.vcf snvs_{pid}_somatic_snvs_conf_8_to_10.VEP.vcf

vep_species="homo_sapiens"
vep_assembly="GRCh37"
vep_out_format="vcf"

input_vcf=${1}
output_vcf=${2}
threads=${VEP_FORKS}

## Annotate the high confidence variants
## Parse for the functional consequences
${PERL_BINARY} ${VEP_SW_PATH} \
  --input_file $input_vcf \
  --species $vep_species \
  --assembly $vep_assembly \
  --output_file STDOUT \
  --format $vep_out_format \
  --fork $threads \
  --fasta ${VEP_FA_INDEX} \
  --everything \
  --vcf \
  --cache \
  --offline \
  --force_overwrite \
  --no_stats \
  --dir_cache ${VEP_CACHE_BASE} \
  | ${PYTHON_BINARY} ${TOOL_PARSE_VEP} | ${BGZIP_BINARY} -f > $output_vcf