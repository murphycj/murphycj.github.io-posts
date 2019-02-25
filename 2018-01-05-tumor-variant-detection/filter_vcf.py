"""
Given a set of parameters for identifying somatic mutations, this script
will
"""

import argparse

import vcf
import pandas

CONTROLS = [
    'HL-WES-12', 'HL-WES-13', 'HL-WES-14', 'HL-WES-15', 'HL-WES-16',
    'HL-WES-39', 'HL-WES-40', 'HL-WES-41', 'HL115-2374', 'HL116-2453',
    'HL117-2598'
]


def main(args):

    support = pandas.read_csv(args.support, index_col=0)
    coverage = pandas.read_csv(args.coverage, index_col=0)

    samples = pandas.read_excel(args.samples, 'samples')

    primary_tumors = list(set(samples['mouse'].tolist()))
    sample_map = dict(zip(
        samples['sample_ID'].tolist(), samples['mouse'].tolist()))

    total_mutations = 0
    proportion_in_dbsnp = 0.0
    proportion_in_dbsnp_2 = 0.0
    proportion_high_quality_singleton = 0.0
    proportion_10 = 0.0
    proportion_50 = 0.0
    proportion_100 = 0.0

    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    count = 0

    for v in vcf_in:
        count += 1
        if count % 100 == 0:
            print(count)

        name = str(v.CHROM) + '-' + str(v.POS)

        # filter by control

        num_with_min_coverage = 0
        filter_because_control = False

        for s in CONTROLS:
            if coverage.ix[name][s] >= args.Cc:
                num_with_min_coverage += 1

            if support.ix[name][s] >= args.Cr:
                filter_because_control = True

        if num_with_min_coverage < args.Cn:
            filter_because_control = True

        # filter by tumor

        filter_because_tumor = False

        primary_tumor_count = {i: False for i in primary_tumors}
        mutant_samples = []

        for s in v.samples:

            if s.sample not in CONTROLS:

                if s.called:
                    primary_tumor_count[sample_map[s.sample]] = True
                    mutant_samples.append(s.data.AD)

                if support.ix[name][s.sample] >= args.Tr:
                    primary_tumor_count[sample_map[s.sample]] = True

        n_primary = sum(map(
            lambda x: 1 if x else 0, primary_tumor_count.values()))

        if n_primary >= args.Tn:
            filter_because_tumor = True

        if not filter_because_tumor and not filter_because_control:
            total_mutations += 1

            if len(mutant_samples) == 1 and mutant_samples[0] >= 20:
                proportion_high_quality_singleton += 1

            if v.ID is not None:
                proportion_in_dbsnp += 1

            if v.ID is not None and n_primary >= 2:
                proportion_in_dbsnp_2 += 1

    fout = open(args.out, 'w')
    fout.write(
        'total_mutation,proportion_in_dbsnp,proportion_in_dbsnp_2,' +
        'proportion_high_quality_singleton,proportion_10,proportion_50,' +
        'proportion_100\n'
    )
    fout.write(str(total_mutations))
    fout.write(',' + str(proportion_in_dbsnp/float(total_mutations)))
    fout.write(',' + str(proportion_in_dbsnp_2/float(total_mutations)))
    fout.write(',' + str(
        proportion_high_quality_singleton/float(total_mutations)))
    fout.write(',' + str(proportion_10/float(total_mutations)))
    fout.write(',' + str(proportion_50/float(total_mutations)))
    fout.write(',' + str(proportion_100/float(total_mutations)))
    fout.close()


parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf', type=str, help='VCF file', required=True)
parser.add_argument('--samples', type=str, help='samples', required=True)
parser.add_argument('--Tn', type=int, help='', required=True)
parser.add_argument('--Tr', type=int, help='', required=True)
parser.add_argument('--Cn', type=int, help='', required=True)
parser.add_argument('--Cc', type=int, help='', required=True)
parser.add_argument('--Cr', type=int, help='', required=True)
parser.add_argument('--coverage', type=str, help='', required=True)
parser.add_argument('--support', type=str, help='', required=True)
parser.add_argument('--out', type=str, help='out results', required=False)
args = parser.parse_args()

main(args=args)
