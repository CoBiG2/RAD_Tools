#!/bin/bash

set -e

vcffile=${1}

function get_individuals() {
  local input=${1}
  _a=( $( head -n12 ${input} |tail -n 1) )
  individuals="${_a[@]:9}"
  echo ${individuals}
  }


function get_counts() {
  local vcf_file="${1}"
  shift 1
  local individuals=( $@ )
  mkdir -p indiv_vcfs
  mkdir -p indiv_counts
  for indiv in "${individuals[@]}"
    do
    vcftools --vcf ${vcf_file} --indv ${indiv} --out indiv_vcfs/${indiv} --recode --max-alleles 2
    vcftools --vcf indiv_vcfs/${indiv}.recode.vcf --counts2 --out indiv_counts/${indiv} 
    done
  }

function transpose_and_join() {
  local vcf_file=${1}
  local countfile=${vcf_file}.counts
  local sizefile=${vcf_file}.size
  shift 1
  local individuals=( $@ )
  local loci=$(tail -n +2 indiv_counts/${individuals[0]}.frq.count| cut -f 1,2 |tr "\t" "."|tr "\n" "\t")
  echo ${loci} > ${countfile} 
  echo ${loci} > ${sizefile}
  for indiv in "${individuals[@]}"
    do
    count=$(tail -n +2 indiv_counts/${indiv}.frq.count|cut -f 5|tr "\n" "\t")
    name_count="${indiv} ${count}"
    size=$(tail -n +2 indiv_counts/${indiv}.frq.count|cut -f 4|tr "\n" "\t")
    name_size="${indiv} ${size}"
    echo ${name_count} >> ${countfile}
    echo ${name_size} >> ${sizefile}
    done
  }

indivs=$(get_individuals ${vcffile})
get_counts ${vcffile} ${indivs}
transpose_and_join ${vcffile} ${indivs}
