import re
import sys
import pandas
import numpy as np
import os
import argparse

def main(args):

    results = []

    for f in args.files:

        temp = os.path.split(f)[1].replace('.csv','').split('_')
        temp = map(lambda x: int(x.split('-')[1]),temp)

        temp = temp + pandas.read_csv(f).ix[0].tolist()

        results.append(temp)

    data = pandas.DataFrame(
        results,
        columns=[
            'Tn','Tr','Cn','Cc','Cr',
            'total_mutations','proportion_in_dbsnp','proportion_in_dbsnp_2',
            'proportion_high_quality_singleton','proportion_10',
            'proportion_50','proportion_100'
        ]
    )

    data.to_csv(args.out,index=False)


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--files',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of cufflinks output files.'
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)
