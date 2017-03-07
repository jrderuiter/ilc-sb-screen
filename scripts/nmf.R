suppressMessages({
    library(argparser)
    library(NMF)
})


# Parse arguments.
parser <- arg_parser('A text file modifying program')

parser = add_argument(parser, '--counts', help='Raw expression counts')
parser = add_argument(parser, '--norm_counts',
                      help='Normalized expression counts')
parser = add_argument(parser, '--output_dir', help='Output directory')
parser = add_argument(parser, '--cores', help='Number of cores to use',
                      type='integer', default=1)

argv = parse_args(p)

# Check required arguments.
if (all(is.na(argv$counts))) {
    stop('Missing required argument --counts')
} else if (is.na(argv$norm_counts)) {
    stop('Missing required argument --norm_counts')
} else if (is.na(argv$output_dir)) {
    stop('Missing required argument --output_dir')
}

# Create directory if it doesn't exist.
if (dir.exists(argv$output_dir)) {
    stop('Output directory already exists')
} else {
    dir.create(argv$output_dir, recursive=TRUE)
}

# Read input counts/expr.
counts = read.delim(argv$counts, sep='\t',
                    check.names=FALSE, row.names=1)
norm_counts = read.delim(argv$norm_counts, sep='\t',
                         check.names=FALSE, row.names=1)

# Mask any genes with too little counts or no variance.
norm_counts = norm_counts[(apply(counts >= 10, 1, sum) >= 2) &
                          (apply(counts, 1, var) > 1e-10), ]

# Base all counts at zero.
norm_counts = norm_counts - apply(norm_counts, 1, min)

# Perform rank estimation.
estimated_rank = nmf(as.matrix(norm_counts), rank=2:5,
                     nrun=30, .pbackend=n_cores)
saveRDS(estimated_rank, file=file.path(argv$output_dir, 'rank_estimation.rds'))

# Do the actual factorization with 4 clusters.
nmf4 = nmf(norm_counts, rank=4, nrun=200, .pbackend=n_cores)
saveRDS(nmf4, file=file.path(argv$output_dir, 'factorization.rds'))

# Extract results of factorization.
nmf4 = readRDS(factorization_path)

nmf4_coef = t(as.data.frame(coef(nmf4)))
colnames(nmf4_coef) = c(1, 2, 3, 4)

write.table(nmf4_coef, file=file.path(argv$output_dir, 'coefficients.txt'),
            sep='\t', row.names=TRUE, quote=FALSE)

# Extract features/genes.
nmf4_features = extractFeatures(nmf4)
nmf4_genes = do.call('rbind',
    lapply(seq_along(nmf4_features), function(i) {
        gene_symbols = rownames(nmf4)[nmf4_features[[i]]]
        data.frame(cluster=i, gene_symbol=gene_symbols)
    })
)

write.table(nmf4_genes, file=file.path(argv$output_dir, 'genes.txt'),
            sep='\t', row.names=FALSE, quote=FALSE)
