"""Script that extracts expression counts from Firebrowse expression file."""

import argparse

import pandas as pd


def main():
    """Main function."""

    args = parse_args()

    # Read counts and swap header levels, so that 'normalized_count' etc.
    # are now top-level instead of the sample names.
    data = pd.read_csv(args.input, sep='\t', header=[0, 1], index_col=0)
    data = data.swaplevel(0, 1, axis=1)

    # Extract the required counts.
    data = data['normalized_count']

    # Extract gene symbols if needed.
    if args.symbols:
        data.index = [i.split('|')[0] for i in data.index]
        data = data.drop('?')

    data.to_csv(args.output, sep='\t', index=True)


def parse_args():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument('input')
    parser.add_argument('output')

    parser.add_argument('--type', default='normalized_count')
    parser.add_argument('--symbols', default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
