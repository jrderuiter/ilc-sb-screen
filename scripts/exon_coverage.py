#!/usr/bin/env python

import argparse
from os import path

import numpy as np
import pandas as pd
import pysam

from tqdm import tqdm


def main():
    args = parse_args()

    agg_func = np.median if args.median else np.mean

    result = coverage(
        args.bams,
        args.gtf,
        transcript_ids=args.transcript_ids,
        agg_func=agg_func)

    result.to_csv(args.output, sep='\t', index=True)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--bams', nargs='+', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--output', required=True)

    parser.add_argument(
        '--transcript_ids', required=False, default=None, nargs='+')
    parser.add_argument('--median', default=False, action='store_true')
    # parser.add_argument('--verbose', default=False, action='store_true')

    return parser.parse_args()


def coverage(bam_paths,
             gtf_path,
             transcript_ids=None,
             verbose=False,
             agg_func=None):

    # Setup record iterator from gtf file.
    gtf_file = pysam.Tabixfile(gtf_path, parser=pysam.asGTF())
    gtf_records = (rec for rec in gtf_file.fetch() if rec.feature == 'exon')

    if transcript_ids is not None:
        transcript_ids = set(transcript_ids)
        gtf_records = (rec for rec in gtf_records
                       if rec['transcript_id'] in transcript_ids)

    if verbose:
        gtf_records = tqdm(gtf_records, leave=False)

    # Build frame.
    rows = _coverage_gen(bam_paths, gtf_records, agg_func=agg_func)
    index_names = ['transcript_id', 'chr', 'start', 'end', 'strand']

    result = pd.DataFrame.from_records(
        rows, columns=index_names + list(bam_paths))
    result = result.set_index(index_names)

    return result


def _coverage_gen(bam_paths, records, agg_func=None):
    bam_files = [pysam.AlignmentFile(fp) for fp in bam_paths]

    for rec in records:
        index = (rec['transcript_id'], rec.contig, rec.start, rec.end,
                 rec.strand)
        region = index[1:-1]

        coverages = tuple(region_coverage(bam_file, region, agg_func=agg_func)
                          for bam_file in bam_files)  # yapf: disable

        yield index + coverages


def region_coverage(bam_file, region, agg_func=None):
    agg_func = agg_func or np.mean
    contig, start, end = region

    # Get coverage within region.
    pileups = bam_file.pileup(contig, start, end, stepper='all', truncate=True)
    coverage = [p.n for p in pileups]

    # Handle empty case.
    if len(coverage) == 0:
        coverage = [0]

    return agg_func(coverage)


if __name__ == '__main__':
    main()
