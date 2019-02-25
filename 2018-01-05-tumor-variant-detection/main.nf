#!/usr/bin/env nextflow

// Pipeline for parameterizing my somatic mutation calling.
// I don't include all data files needed to reproduce this pipeline, but if
// you really want the data then please ask.

// parameters and tools for pipeline

basedir = "$baseDir"

params.bgzip = '/home/chm2059/chm2059/lib/tabix-0.2.6/bgzip'
params.tabix = '/home/chm2059/chm2059/lib/tabix-0.2.6/tabix'
params.bcftools = '~/chm2059/lib/bcftools-1.3.1/bcftools'
params.vcftools = '~/chm2059/lib/vcftools_0.1.12b/bin/vcftools'
params.vcftools_bin = '~/chm2059/lib/vcftools_0.1.12b/bin/'
params.vcf2bed = '~/chm2059/lib/bedops_linux_x86_64-v2.4.2/vcf2bed'
params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'
params.varscan = '/home/chm2059/chm2059/lib/VarScan.v2.3.7.jar'
params.gatk = '~/chm2059/lib/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
params.reference = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.reference_dict = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.reference_fai = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.dbsnp_indels = '/home/chm2059/chm2059/data/refdata/mm10/mgp.v4.indels.dbSNP.vcf'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.gtf'
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.snpsift = '/home/chm2059/chm2059/lib/snpEff/SnpSift.jar'
params.snpeff = '/home/chm2059/chm2059/lib/snpEff/snpEff.jar'
params.samtools = '~/chm2059/lib/samtools-1.2/samtools'

params.combinepileup = 'collect_stats.py'
params.filtervcf = 'filter_vcf.py'
params.aggregate = 'aggregate_results.py'
params.samples = 'samples.xlsx'

// contains list of somatic variants called on all tumors using some very
// non-stringent parameters

initial_somatic = file('somatic.varscan.vcf')

// parameters to test over

Tn = Channel.from(2,3,4,30)
Tr = Channel.from(2,3,4)
Cn = Channel.from(1,2,4,11)
Cc = Channel.from(3,5,10)
Cr = Channel.from(2,3,4)

// annotation somatic mutations with DBSNP

process dbsnp {

  storeDir "${baseDir}"

  input:
    file initial_somatic

  output:
    set file('somatic.varscan.dbsnp.vcf') into dbsnp_out_vcf2bed
    set file('somatic.varscan.dbsnp.vcf') into dbsnp_out_filter
    set file('somatic.varscan.dbsnp.vcf') into dbsnp_out_pileup

  """
  ${params.java} -Xmx12g -jar \
    ${params.snpsift} annotate \
    ${params.dbsnp} \
    ${initial_somatic} \
    | ${params.java} -Xmx16g -jar \
    ${params.snpsift} annotate \
    ${params.dbsnp_indels} \
    | ${params.java} -Xmx16g -jar \
    ${params.snpeff} eff \
    -no-downstream \
    -no-upstream \
    -v GRCm38.74 > somatic.varscan.dbsnp.vcf
  """
}

// convert annotated vcf file to bed file

process vcf_to_bed {

  storeDir "${baseDir}"

  input:
    file dbsnp_out_vcf2bed

  output:
    set file('variant_positions.bed') into variant_positions_wex
    set file('variant_positions.bed') into variant_positions_rnaseq

  """
  export PATH=~/chm2059/lib/bedops_linux_x86_64-v2.4.2/:$PATH
  ${params.vcf2bed} \
     < ${dbsnp_out_vcf2bed} \
    | awk -v OFS='\t' '{print \$1, \$2, \$3}' \
    | uniq > \
    variant_positions.bed
  """
}

// generate the list of all normal and tumor NGS samples I have

initial_rnaseq = (01..114).collect{
  if (it < 10) {
    'HL0' + it
  } else {
    'HL' + it
  }
}
initial_rnaseq << "HL115-2374"
initial_rnaseq << "HL116-2453"
initial_rnaseq << "HL117-2598"

initial_wex = (01..16).collect{
  if (it < 10) {
    'HL-WES-0' + it
  } else {
    'HL-WES-' + it
  }
}
initial_wex << "HL-WES-39"
initial_wex << "HL-WES-40"
initial_wex << "HL-WES-41"


// get the list of BAM/BAI files for all tumor and normal samples

rnaseq_bams = Channel.empty()
wex_bams = Channel.empty()
bams = Channel.empty()
test = Channel.empty()

Channel
  .fromFilePairs("/home/chm2059/chm2059/elementoCantley_2014_9_29/results/RNAseq/STARGATK/HL*/HL*.gatk.bam", size: 1)
  .filter{
    it[0] in initial_rnaseq
  }
  .flatMap{
    it -> [[it[0], it[1][0]]]
  }
  .spread(variant_positions_rnaseq)
  .set{
    rnaseq_bams
  }

Channel
  .fromFilePairs("/home/chm2059/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/HL*/HL*.gatk.bam", size: 1)
  .filter{
    it[0] in initial_wex
  }
  .flatMap{
    it -> [[it[0], it[1][0]]]
  }
  .spread(variant_positions_wex)
  .set{
    wex_bams
  }

test.concat( rnaseq_bams, wex_bams ).set{bams}

// generate pileups for all tumor and normal samples over the locations
// of the identifed somatic mutations

process get_pileup {

  storeDir "${baseDir}/mpileups/${s}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=4G -pe smp 2 -l h_rt=1:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set s, bam, file(variant_positions_file) from bams

  output:
    set s into all_pileup_samples
    set file('*.pileup') into all_pileup_files

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  rsync -L ${bam} tmp.bam
  rsync -L ${bam}.bai tmp.bam.bai
  ${params.samtools} mpileup \
    -f ${params.reference} \
    -l ${variant_positions_file} \
    -Q 15 \
    tmp.bam \
    -o ${s}.pileup
  """
}

// combine the pileups into two files: one containing the number of reads
// that support the somatic mutation and another that contains the
// total number of reads overlapping each variant position

process combine_pileup {

  storeDir "${baseDir}"

  input:
    file ap from all_pileup_files.toList()
    val s from all_pileup_samples.toList()
    file dbsnp_out_pileup

  output:
    set file('coverage.csv') into combine_pileup_coverage
    set file('support.csv') into combine_pileup_support

  """
  python ${params.combinepileup} \
    --vcf ${dbsnp_out_pileup} \
    --pileups ${ap} \
    --samples ${s.join(' ')}
  """
}

// using the vcf of somatic variants and combined pileup files, iterate
// over all parameter combinations to compute various metrics of the
// algorithm's ability to identify true somatic variants

filter_in = dbsnp_out_filter
  .spread(combine_pileup_coverage)
  .spread(combine_pileup_support)
  .spread(Tn)
  .spread(Tr)
  .spread(Cn)
  .spread(Cc)
  .spread(Cr)

process filter {

  storeDir "${baseDir}/filterStats/Tn-${Tn}_Tr-${Tr}_Cn-${Cn}_Cc-${Cc}_Cr-${Cr}"
  maxForks 20

  input:
    set file(vcf), file(coverage), file(support), Tn, Tr, Cn, Cc, Cr from filter_in

  output:
    set file('*.csv') into filter_out

  """
  python ${params.filtervcf} \
    --vcf ${vcf} \
    --samples ${params.samples} \
    --Tn ${Tn} \
    --Tr ${Tr} \
    --Cn ${Cn} \
    --Cc ${Cc} \
    --Cr ${Cr} \
    --coverage ${coverage} \
    --support ${support} \
    --out Tn-${Tn}_Tr-${Tr}_Cn-${Cn}_Cc-${Cc}_Cr-${Cr}.csv

  """
}

// combine the result

process combine_results {

  storeDir "${baseDir}"

  input:
    file combine_files from filter_out.toList()

  output:
    set file('results.csv') into combine_results_out

  """
  python ${params.aggregate} --files ${combine_files} --out results.csv
  """
}
