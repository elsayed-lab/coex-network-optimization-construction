#!/bin/env Rscript
#
# Co-expression Network Parameter Optimization (v6)
# Keith Hughitt (khughitt@umd.edu)
# 2017/04/06
#
# Usage: coex_network_param_opt.R <network_type> <low_count_threshold> 
#                                 <cpm> <log2> <voom>
#                                 <quantile_normalize> <batch_adjustment>
#                                 <similarity_measure> <adj_power> 
#                                 <topological_overlap> <merge_cor>
#
# All of the parameters must be included and must be specified in the above
# ordering.
#
# Parameters:
#
#   network_type        [signed|unsigned]
#   low_count_threshold int
#   cpm                 boolean
#   log2                boolean
#   voom                boolean
#   quantile_normalize  boolean
#   batch_adjustment    [limma|combat|none]
#   similarity_measure  [cor|bicor|spearman|dist]
#   adj_power           int
#   topological_overlap boolean
#   merge_cor           0-1
#
# For a more detailed explanation of each parameter, see the main Co-expression
# Network analysis.


# Set working directory manually to be safe
setwd(file.path(Sys.getenv('RESEARCH'), '2015/115-coex-network-param-opt-v6'))

library(knitr)
knit('../00-shared/Rmd/init/header_network.Rmd', quiet=TRUE, output=tempfile())

sessionInfo()

opts_knit$set(progress=TRUE, verbose=TRUE)
opts_chunk$set(error=FALSE)

# Step 1: load parameter arguments from command line
args <- commandArgs(trailingOnly=TRUE)

outfile             <- args[1]
settings            <- args[2]
network_type        <- args[3]
low_count_threshold <- as.numeric(args[4])
cpm_transform       <- as.logical(args[5])
log2_transform      <- as.logical(args[6])
voom_transform      <- as.logical(args[7])
quantile_normalize  <- as.logical(args[8])
batch_adjustment    <- args[9]
similarity_measure  <- args[10]
adj_power           <- as.numeric(args[11])
topological_overlap <- as.logical(args[12])
merge_cor           <- as.numeric(args[13])

# Step 2: Read in base settings
knit(settings, output=tempfile())

# Step 3: Overide base settings
CONFIG$network_type               <- network_type
CONFIG$low_count_threshold        <- low_count_threshold
CONFIG$network_cpm                <- cpm_transform
CONFIG$network_log2               <- log2_transform
CONFIG$network_voom               <- voom_transform
CONFIG$network_quantile_normalize <- quantile_normalize
CONFIG$network_batch_adjust       <- batch_adjustment
CONFIG$similarity_measure         <- similarity_measure
CONFIG$adjmatrix_power            <- adj_power
CONFIG$topological_overlap        <- topological_overlap
CONFIG$merge_correlation          <- merge_cor

# Disable resampling
CONFIG$enable_resampling <- FALSE

CONFIG$debug         <- TRUE
CONFIG$verbose       <- TRUE
CONFIG$include_plots <- FALSE
CONFIG$use_cache     <- FALSE

# Load defaults
knit('../00-shared/Rmd/init/load_settings_network.Rmd', quiet=TRUE, output=tempfile())

print("===========================================================")
print("Co-expression Network Construction Settings:")
print(sprintf("%20s: %s", 'Settings file', settings))
print(sprintf("%20s: %s", 'Network type', network_type))
print(sprintf("%20s: %s", 'Low count threshold', low_count_threshold))
print(sprintf("%20s: %s", 'CPM transform', cpm_transform))
print(sprintf("%20s: %s", 'Log2 transform', log2_transform))
print(sprintf("%20s: %s", 'Voom transform', voom_transform))
print(sprintf("%20s: %s", 'Quantile normalize', quantile_normalize))
print(sprintf("%20s: %s", 'Batch adjustment', batch_adjustment))
print(sprintf("%20s: %s", 'Similarity measure', similarity_measure))
print(sprintf("%20s: %s", 'Adjacency power', adj_power))
print(sprintf("%20s: %s", 'Topological overlap', topological_overlap))
print(sprintf("%20s: %s", 'Module merge correlation', merge_cor))
print("===========================================================")

# Print variable creations command to make it easier to reproduce
print(sprintf("args <- c('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')",
              outfile, settings, network_type, low_count_threshold, cpm_transform,
              log2_transform, voom_transform, quantile_normalize, batch_adjustment,
              similarity_measure, adj_power, topological_overlap, merge_cor))

# Main analysis
knit('../00-shared/Rmd/init/load_counts.Rmd', quiet=TRUE, output=tempfile())

if (CONFIG$target == 'host') {
    if (CONFIG$host == 'D. melanogaster') {
        knit('../00-shared/Rmd/init/load_fly_annotations.Rmd', quiet=TRUE, output=tempfile())
    } else if (CONFIG$host == 'C. elegans') {
        knit('../00-shared/Rmd/init/load_worm_annotations.Rmd', quiet=TRUE, output=tempfile())
    } else {
        knit('../00-shared/Rmd/init/load_host_annotations.Rmd', quiet=TRUE, output=tempfile())
    }
} else {
    knit('../00-shared/Rmd/init/load_pathogen_annotations.Rmd', quiet=TRUE, output=tempfile())
}
knit('../00-shared/Rmd/init/create_expression_set.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/main/filter_counts.Rmd', quiet=TRUE, output=tempfile())

if (CONFIG$target == 'pathogen' && CONFIG$filter_multicopy) {
    knit('../00-shared/Rmd/main/filter_multicopy_genes.Rmd', quiet=TRUE, output=tempfile())
}

knit('../00-shared/Rmd/main/data_prep_de.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/main/differential_expression.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/main/data_prep_network.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/main/network_construction.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/main/network_module_detection.Rmd', quiet=TRUE, output=tempfile())


knit('../00-shared/Rmd/results/overview.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/results/expression_profiles.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/results/go_enrichment_network.Rmd', quiet=TRUE, output=tempfile())
knit('../00-shared/Rmd/results/kegg_enrichment_network.Rmd', quiet=TRUE, output=tempfile())

# Host-specific results
if (CONFIG$target == 'host' && CONFIG$host %in% c('H. sapiens', 'M. musculus')) {
    # iRefIndex / ConsensusPathDB
    knit('../00-shared/Rmd/results/host.Rmd', quiet=TRUE, output=tempfile())

    if (CONFIG$host == 'H. sapiens') {
        # Marbach et al. (2016) TF target enrichment
        knit('../00-shared/Rmd/results/hsapiens_marbach2016_tf_enrichment.Rmd', quiet=TRUE, output=tempfile())
    }
}

# Summarize module annotation enrichment
go_summary <- summarize_enrichment_result(module_go_enrichment)
kegg_summary <- summarize_enrichment_result(module_kegg_enrichment)

# total of differentially expressed genes
total_de <- sum(num_de)

# median module size
med_module_size <- as.numeric(median(table(module_labels)))

# total number of genes after filtering
num_genes <- length(gene_ids)

# param opt id
test_id <- sub('.out', '', basename(outfile))

# create result vector
entries <- c(test_id, network_type, low_count_threshold, cpm_transform,
             log2_transform, voom_transform, quantile_normalize,
             batch_adjustment, similarity_measure, adj_power,
             topological_overlap, merge_cor, num_modules_before, num_modules,
             num_genes, med_module_size, go_summary$total_categories,
             go_summary$unique_categories, go_summary$num_enriched_modules,
             go_summary$mean_pval, kegg_summary$total_categories,
             kegg_summary$unique_categories, kegg_summary$num_enriched_modules,
             kegg_summary$mean_pval, total_de)

# For human / mouse, also include enrichment and overlap with iRefIndex, CPDB,
# and Marbach et al. (2016) networks
if (CONFIG$target == 'host' && CONFIG$host %in% c('H. sapiens', 'M. musculus')) {
    # Overlap with iRefIndex PPI network
    total_edge_weight     <- sprintf("%0.2f", sum(adjacency_matrix))
    irefindex_edge_weight <- sprintf("%0.2f", iref_edge_sum)
    entries <- append(entries, c(total_edge_weight, irefindex_edge_weight))
    
    # CPDB enrichment
    cpdb_summary <- summarize_enrichment_result(cpdb_pathway_enrichment)
    entries <- append(entries,
                     c(cpdb_summary$total_categories, cpdb_summary$unique_categories,
                       cpdb_summary$num_enriched_modules, cpdb_summary$mean_pval))

    if (CONFIG$host == 'H. sapiens') {
        # Marbach et al. (2016) TF co-regulated enrichment (Human only)
        marbach_summary <- summarize_enrichment_result(module_coreg_enrichment)
        entries <- append(entries,
                         c(marbach_summary$total_categories, marbach_summary$unique_categories,
                           marbach_summary$num_enriched_modules, marbach_summary$mean_pval))
    }
}

# Output results table entry
output_row <- paste0(paste(entries, collapse=','), '\n')
cat(output_row, file=outfile, sep='')

# Save gene module assignments
module_assignments_file <- sub('\\.out', '.RData', outfile)

# Save module mapping
module_mapping <- result %>% select(gene_id, color)
save(module_mapping, file=module_assignments_file) 

# Save rounded adjacency matrix
# 2016/09/12 - disabling until explicitly needed (requires significant storage)
#adjacency_matrix_file <- sub('\\.out', '_adjmat.RData', outfile)
#adjacency_matrix <- round(adjacency_matrix[upper.tri(adjacency_matrix)], 3)
#save(adjacency_matrix, file=adjacency_matrix_file)

print("Results:")
print(sprintf("num modules: %d", num_modules))
print(sprintf("num genes: %d", num_genes))
print(sprintf("num enriched GO: %d", go_summary$total_categories))
print(sprintf("num enriched KEGG: %d", kegg_summary$total_categories))
date()

