"""
filter_tumor_from_control_set.py

Given a VCF file containing tumor and normal (control) samples and their
corresponding set of BAM files, remove mutations from the "somatic" call set
using the given parameters. Tailored for RNA-seq data

"""

import vcf
import multiprocessing
import sys
import pandas
import os
import argparse


class Pileup:
    """
    store data on genomic variations (snps and small indels)
    """

    def __init__(self, filename, pileup_data):

        self.seq_id = pileup_data[0]
        self.position = int(pileup_data[1])
        self.reference_base = pileup_data[2]
        self.coverage = int(pileup_data[3])

        #read bases and qualities are not printed if coverage is 0

        if self.coverage==0:
            self.read_bases = ''
            self.read_qualities = ''
        else:
            if len(pileup_data)!=6:
                print pileup_data
                print filename
                sys.exit('Did not find 6 columns in pileup data!')
            self.read_bases = pileup_data[4]
            self.read_qualities = pileup_data[5]

    def base_count(self, base):
        count = self.read_bases.count(base.lower()) + self.read_bases.count(base.upper())
        if (base.upper() == self.reference_base) or (base.lower() == self.reference_base):
            count += self.read_bases.count('.')
            count += self.read_bases.count(',')
        return count

def parse_vcf(vcf_file):

    mutant_sites = {}

    vcf_in = vcf.Reader(open(vcf_file,'r'))
    for v in vcf_in:

        alternative_alleles = []

        if len(v.REF)==1:

            for alternative in v.ALT:
                if len(alternative)>1:
                    alternative = str(alternative)[1::]

                    alternative = '+' + str(len(alternative)) + alternative

                    alternative_alleles.append(alternative)
                else:
                    alternative_alleles.append(str(alternative))
        else:
            reference = str(v.REF)[1::]

            reference = '-' + str(len(reference)) + reference

            alternative_alleles.append(reference)

        mutant_sites[str(v.CHROM) + '-' + str(v.POS)] = alternative_alleles
    return mutant_sites

def get_filter_site_set(args,pileup_files,mutant_sites,coverage,support):

    for pileup_file, sample in pileup_files:
        for line in open(pileup_file,'r').readlines():
            pileup = Pileup(pileup_file,line.rstrip().split('\t'))

            site = str(pileup.seq_id) + '-' + str(pileup.position)


            coverage.ix[site][sample]=pileup.coverage
            mutant_support_count=0
            for b in mutant_sites[site]:
                mutant_support_count += pileup.base_count(b)

            support.ix[site][sample]=mutant_support_count


    support = support.fillna(0)
    coverage = coverage.fillna(0)

    coverage.to_csv('coverage.csv')
    support.to_csv('support.csv')

def main(args):

    mutant_sites = parse_vcf(args.vcf)

    coverage = pandas.DataFrame(index=mutant_sites.keys(),columns=args.samples)
    support = pandas.DataFrame(index=mutant_sites.keys(),columns=args.samples)

    pileup_files = zip(args.pileups,args.samples)

    get_filter_site_set(args,pileup_files,mutant_sites,coverage,support)


parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,required=True)
parser.add_argument('--pileups',type=str,required=True,help='Pileup files for control samples',nargs='+')
parser.add_argument('--samples',type=str,required=True,help='Sample names (same order as pileups)',nargs='+')
args = parser.parse_args()

main(args=args)
