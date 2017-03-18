"""Script that uses scikit-bio to convert fasta files to fastq files."""

import argparse

from skbio import io as skbio_io


def main():
    """Main function."""

    args = parse_args()

    if args.output.endswith('.gz'):
        compression = 'gzip'
    else:
        compression = 'auto'

    with skbio_io.open(args.output, mode='w', compression=compression) as file_:
        seqs = skbio_io.read(args.fasta, format='fasta', qual=args.qual)

        if args.qual is None:
            for seq in seqs:
                # Missing quality, add default 'dummy' quality scores.
                seq.positional_metadata['quality'] = args.default_qual
                skbio_io.write(seq, into=file_, format='fastq',
                               variant=args.variant)
        else:
            for seq in seqs:
                skbio_io.write(seq, into=file_, format='fastq',
                               variant=args.variant)


def parse_args():
    """Parses arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument('fasta')
    parser.add_argument('output')

    parser.add_argument('--qual', default=None)
    parser.add_argument('--default_qual', type=int, default=40)
    parser.add_argument('--variant', default='sanger')

    return parser.parse_args()


if __name__ == '__main__':
    main()
