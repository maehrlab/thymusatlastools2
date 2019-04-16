## ------------------------------------------------------------------------

# Several different methods of variable gene selection, all based on outliers in mean-versus-sd plots.
# You can give either x and y cutoffs for these plots or the number of genes to select -- not both.
#' @export
var_gene_select = function( dge, results_path, 
                            excess_var_cutoff = 0.5,
                            log_expr_cutoff = 0.1 ){
  
  # # Need to do this to initialize @mean.var even if not going to use those cutoffs
  dge = Seurat::FindVariableGenes( dge, x.low.cutoff = log_expr_cutoff, y.cutoff = excess_var_cutoff)

  # # Save variable genes and parameters
  vgsrp = file.path( results_path, "var_gene_select" )
  dir.create.nice( vgsrp )
  cell_markers = c()
  variable_cell_markers = intersect( cell_markers, as.character( dge@var.genes ) )
  variable_cell_markers = c(paste0(length(variable_cell_markers), "total"), variable_cell_markers)
  text2file(file.path(vgsrp, "markers_among_variable_genes.txt"), variable_cell_markers)
  totalstring = paste(length(as.character(dge@var.genes)), "total")
  var_genes = c(totalstring, as.character(dge@var.genes))
  text2file(file.path(vgsrp, "variable_genes.txt"), var_genes)
  if( is.null( excess_var_cutoff    ) ) { excess_var_cutoff    = "NULL" }
  if( is.null( log_expr_cutoff      ) ) { log_expr_cutoff      = "NULL" }
  vsp        = c(  excess_var_cutoff,   log_expr_cutoff)
  names(vsp) = c( "excess_var_cutoff", "log_expr_cutoff")
  text2file(file.path(vgsrp, "var_gene_selection_params.txt"), collapse_by_name(vsp))
  
  return(dge)
}

#' @export
save_complexity_plot = function(dge, result_path){
  dir.create.nice( file.path( result_path, "QC" ) )
  f = file.path(result_path, "QC/complexity.pdf")
  p = ggplot( data.frame( complexity = colSums( dge@raw.data > 0 ) ) ) + 
    ggtitle("Complexity") + geom_histogram(aes(x = complexity)) 
  ggsave(f, p)
  f = file.path(result_path, "QC/UMIs_by_cell.pdf")
  p = ggplot( data.frame( UMIs_by_cell = colSums( dge@raw.data ) ) ) + 
    ggtitle("UMIs per cell") + geom_histogram(aes(x = log10(UMIs_by_cell)))
  ggsave(f, p)
}


## ------------------------------------------------------------------------
#' Label doublets using the DoubletFinder package.
#'
#' @param dge seurat object.
#' @param orig.ident Replicate. Only looks for doublets within reps.
#' @param results_path Where to save a brief table summarizing results. If null, nothing gets saved. 
#'
#' @export
LabelDoublets = function( dge, rep_field = "orig.ident", results_path = NULL ){
  data(doublet_rates)
  doublet_rates$rate_pct %<>% as.character %>% gsub("%", "", .) %>% as.numeric
  rate_model = lm( rate_pct ~ recovered_cells, data = doublet_rates )
  reps = dge %>% FetchData(rep_field) %>% extract2(1) %>% unique %>% as.character
  dubs_annot = as.list(reps) %>% setNames( reps )
  for( my_rep in reps ){
    dge_rep = SubsetDataFlex( dge, rep_field, rep_field %>% paste0( " == my_rep" ) )
    recovered_cells = dge_rep@cell.names %>% length
    doublet_rate_pct = predict( rate_model, newdata = data.frame( recovered_cells = recovered_cells ) )
    dge_rep %<>% DoubletFinder::doubletFinder(nExp = recovered_cells * doublet_rate_pct / 100, pK = 0.01 )
    dubs_annot[[my_rep]] = data.frame(pANN = dge_rep@meta.data$pANN_,
                                      is_doublet = as.numeric(dge_rep@meta.data$DF.classifications_=="Doublet"))
    dubs_annot[[my_rep]] %<>% set_rownames( dge_rep@cell.names )
  }
  dubs_annot_catted = Reduce( dubs_annot, f = rbind )
  dge %<>% AddMetaData( dubs_annot_catted )
  counts = dge %>%
    FetchData( c(rep_field, "is_doublet") ) %>%
    table 
  counts %<>% rbind( colSums(counts) ) 
  rownames(counts)[nrow(counts)] = "total"
  cat( "\n\nDoublet counts:\n" )
  print( counts )
  cat( "\n\n" )
  if( !is.null( results_path ) ){
    write.table(counts, file.path(results_path, "doublet_counts.tsv"))
  }
  return( dge )
}

#' Decide if a gene is a ribosomal subunit or a pre-rRNA transcript.
#'
#' @export
#'
is_rp = function( genes ){
  genes %<>% toupper
  p1 = (genes == "RN45S")
  p2 = substring(genes, 1, 3) %in% c("RPS", "RPL")
  return(p1|p2)
}


#' Add up heat shock protein and chaperonin transcripts cell by cell.
#'
#' @export
#'
AddTotalHSP = function( dge ){
  hsp_genes = dge %>% AvailableData %>% grep("^HSP|^CCT", ., ignore.case = T, value = T)
  dge@meta.data$HSP_total = rowSums(FetchData(dge, hsp_genes))
  hsp_genes_plus = c( hsp_genes, "MAGED2", "EPCAM" )
  dge@meta.data$HSP_signature = rowSums(FetchData(dge, hsp_genes_plus))
  return( dge )
}


#' Decide if genes are mitochondrial. Currently uses a crude "grep mt".
#'
#' @export
#'
is_mt = function( genes ){
  genes %<>% toupper
  p = substring(genes, 1, 2) %in% "MT"
  return(p)
}



#' Quantify ribosomal transcript content cell by cell.
#'
#' @export
#'
add_rp_percentage = function(dge) add_rp_mt_percentage(dge )

#' Quantify ribosomal and mitochondrial transcript content cell by cell.
#'
#' @export
#'
add_rp_mt_percentage = function( dge ){
  nUMI = dge %>% FetchData("nUMI") %>% extract2(1) 
  sum_umis = function(predicate){
    raw = dge@raw.data[rownames(dge@data), colnames(dge@data)]
    genes = rownames( raw )
    x = dge@raw.data[ predicate( genes ), ] %>% (Matrix::colSums)
    atat( all( x < nUMI ) ) 
    return(x)
  }
  nUMI_rp = sum_umis(predicate = is_rp)
  nUMI_mt = sum_umis(predicate = is_mt)
  to_add = data.frame( nUMI_mt = nUMI_mt,
                       nUMI_rp = nUMI_rp, 
                       nUMI_pct_mt   = nUMI_mt / nUMI,
                       nUMI_pct_rp   = nUMI_rp / nUMI, 
                       mt_nUMI_pct   = nUMI_mt / nUMI,
                       ribo_nUMI_pct = nUMI_rp / nUMI, 
                       row.names = dge@cell.names )
  dge %<>% AddMetaData( metadata = to_add )
  return(dge)
}


## ------------------------------------------------------------------------

#' Estimate which PC's to keep using the Marchenko-Pastur asymptotic upper bound.
#'
#' @param dge A Seurat object with gene selection and PCA done.
#'
#' @export
#' 
GetDimMarchenkoPastur = function( dge ){
  pc_cutoff = RMTstat::qmp( 1, ndf=length(dge@cell.names), pdim=length(dge@var.genes), var=1)
  is_significant = dge@dr$pca@sdev^2 > pc_cutoff
  num_pc = sum(is_significant)
  return(num_pc)
}


#' Apply a bagged clustering procedure to a Seurat object using PCA as input.
#'
#' @param dge A Seurat object.
#' @param k Number of clusters. (Some may end up empty.)
#' @param num_pc Number of PC's to use from dge. Specify this or feature_names, not both. Default 25. 
#' @param my_seed RNG seed.
#' @param B Number of bootstrap replicates.
#' @param feature_names Fetched via FetchData as input. Specify this or num_pc, not both. Default unused.
#' @param ... Additional args passed to clue::cl_bag.
#'
#' @export
#' 
ClusterBag = function(dge, k, num_pc = NULL, feature_names = NULL, my_seed = 100, B = 100, ... ){
  set.seed(my_seed)
  if(!is.null(feature_names) && !is.null(num_pc) ){
    stop("Please specify only one of 'num_pc' and 'feature_names'.")
  } 
  if(is.null(feature_names) && is.null(num_pc) ){
    warning("Please specify one of 'num_pc' and 'feature_names'. Using 25 PC's by default.")
    num_pc = 25
  }
  if( !is.null( num_pc ) ){
    feature_names = paste0("PC", 1:num_pc)
  } 
  X = FetchData( dge, feature_names )
  
  
  cl_obj = clue::cl_bag( X, k = k, B = B, ... ) 
  memberships = as.data.frame( cl_obj$.Data[,] )
  colnames(memberships) = paste0("cl_bag_", 1:ncol(memberships))
  rownames(memberships) = rownames(dge@meta.data)
  assignment = apply(memberships, 1, which.max)
  confidence = apply(memberships, 1, max)
  memberships$cl_bag_assignment = paste0("cl_bag_", assignment)
  memberships$cl_bag_confidence = confidence
  dge %<>% AddMetaData( memberships )
  return(dge)
}

#' Wrapper for ClusterBag. Measures stability over a range of model sizes.
#'
#' @export
#'
ClusterBagStability = function( dge, k_max = 15, ... ){
  kvals = 2:k_max
  stability = kvals
  for( k in kvals ){
    dge = ClusterBag( dge, k, ...)
    stability[k - 1] = mean(FetchData(dge, "cl_bag_confidence")[[1]])
  }
  plot(kvals, stability)
  return( data.frame(k = kvals, stability = stability) )
}


# # This helps compare two different analyses of the same cells.
# # You put in the usual Seurat object and results path
# # plus another Seurat object that you want to compare against,
# # and names for each of them.
#' @export
compare_views = function(dge, results_path, comparator_dge, dge_name, comparator_name){
  figname = paste0("embedding=", dge_name, "|colors=", comparator_name, ".pdf")
  dge = Seurat::AddMetaData( dge, metadata = comparator_dge@ident[dge@cell.names] , col.name = "other_id" )
  ggsave(file.path(results_path, figname),
         custom_feature_plot(dge, colour = "other_id"),
         width = 5.5, height = 5)
}


#' Rename clusters using an existing metadata field.
#'
#' @param dge : a Seurat object with field `@scale.data` filled in.
#' @param ident.use Cluster labels. The output is the same as this clustering, but renamed. 
#' @param annot.use embryonic day. For example, if you put eday, a cluster might be called "e11_5_c2".
#' @param new.name What field to put the new item in. By default, uses the "ident" slot. 
#'
#' @return character vector.
#' @export
#'
RenameClusters = function( dge, ident.use = "ident", annot.use = "eday", new.name = "ident" ){
  days_present = unique(FetchData(dge, annot.use)[[1]] %>% as.character)
  get_eday_props  = function( x ) x %>% as.character %>% c(days_present) %>% table(exclude=NULL) %>% subtract(1) %>% prop.table
  make_eday_names = function( props ) props %>% is_greater_than(0.2) %>% which %>% names %>% paste0(collapse = "_")
  new_names = aggregate_nice( FetchData( dge, annot.use )[[1]] %>% as.character, 
                              FetchData( dge, ident.use ), 
                              FUN = function(x) x %>% get_eday_props %>% make_eday_names ) 
  old_names = rownames(new_names)
  new_names = c(new_names)
  add_numbers = function( x, sep = "_" ) {
    ave( seq_along(x), x, FUN = function(y) sort(order(y)) ) %>% paste0( x, sep, .  )
  }
  new_names %<>% add_numbers(sep = "_c")
  names(new_names) = old_names
  current_ident = FetchData(dge, ident.use)[[1]] %>% as.character
  if( new.name == "ident"){
    dge %<>% SetIdent(    new_names[current_ident], cells.use = dge@cell.names )
  } else {
    dge %<>% AddMetaData( new_names[current_ident] %>% setNames( dge@cell.names ), col.name = new.name )
  }
  return( dge )
}

## ------------------------------------------------------------------------

#' Find genes with expression patterns similar to the genes you've specified.
#'
#' @param `dge` : a Seurat object with field `@scale.data` filled in.
#' @param `markers`: a character vector; giving gene names.
#' @param `n`: integer; number of results to return.
#' @param `genes.use`: character vector; list of genes eligible to be returned.
#' @param `anticorr` : allow negatively correlated genes; defaults to `FALSE`.
#' @param `aggregator` : If multiple genes, how to combine their correlations. Default is "sum". For an "and"-like filter, use "min".   
#'
#' Given a Seurat object and a list of gene names, this function returns genes 
#' that are strongly correlated with those markers. 
#'
#' @return character vector.
#' @export
#'
get_similar_genes = function( dge, markers, n, genes.use = rownames(dge@data), anticorr = F ){
  # Sanitize input
  if(!all(markers %in% AvailableData(dge))){ 
    warning("Some of your markers have no data available.")
  }
  markers %<>% intersect( AvailableData(dge) ) 
  if(!all(genes.use %in% rownames(dge@data))){ 
    warning("Some of your genes have no data available.")
  }
  genes.use %<>% intersect( rownames(dge@data) ) 
  data.use = dge@data[genes.use, ] %>% Matrix::Matrix( sparse = T )
  assertthat::assert_that(grepl("Matrix", class(data.use)[[1]]))

  # Compute gene-wise mean and sd; optimized for speed
  gene_averages = Matrix::rowMeans( data.use ) 
  N = ncol(data.use)
  squares = Matrix::rowSums(data.use^2)
  gene_vars = squares / N - gene_averages^2
  gene_sds = gene_vars %>% sqrt
  
  # Compute correlations with query
  # Actually, it's halfway between a covariance and a correlation because
  # I divide by the sd of the genes but not the sd of the query
  query = as.matrix( Seurat::FetchData(dge, markers) )
  covariances_all = data.use %*% query - gene_averages %o% colSums(query) 
  correlations_all = Matrix::Diagonal(x = 1 / gene_sds) %*% covariances_all 
  
  # If multiple queries, aggregate by adding. Skim off top of list.
  correlation = Matrix::rowSums( correlations_all )
  correlation = correlation[ setdiff( names( correlation ), markers ) ]
  if( anticorr ){
    similar_genes = names( sort( abs( correlation ), decreasing = T )[ 1:n ] )
  } else {
    similar_genes = names( sort( correlation, decreasing = T )[ 1:n ] )
  }
  return( similar_genes )  
}

# get_similar_genes(thymus_test, "Epcam", 10 )


## ------------------------------------------------------------------------
#' Programmatically feed gene-sets to Enrichr and save formatted results.
#'
#' @export
#'
do_enrichr = function( results_path, geneset, geneset_name, 
                       desired_db = c( "KEGG_2016", 
                                       "WikiPathways_2016",
                                       "Reactome_2016",
                                       "BioCarta_2016",
                                       "Panther_2016",
                                       "NCI-Nature_2016", 
                                       "GO_Biological_Process_2015" ),
                       N_ANNOT_PER_DB = 2 ){
  if(length(geneset) == 0) {
    warning( "List of length zero fed to do_enrichr. Returning NULL.")
    return(NULL)
  }
  
  my_rp = file.path( results_path, paste0("annot_", geneset_name) )
  dir.create.nice( my_rp )
  
  # # Get enrichr results and parse them
  output_table = enrichR::enrichr( genes = geneset, databases = desired_db ) 
  output_table %<>% (base::mapply)(desired_db, FUN = function(X, db) {try({X$database = db}); return(X)}, SIMPLIFY = F)
  output_table %<>% (base::Reduce)( f=rbind )
  output_table %<>% (dplyr::group_by)( database ) %>% (dplyr::top_n)( wt = -P.value, n = N_ANNOT_PER_DB) %>% as.data.frame
  output_table %<>% (dplyr::mutate)( log10_qval = round( log10( Adjusted.P.value ), 1 ) )
  output_table = output_table[c("database", "Term", "Overlap", "log10_qval", "Genes")]
  
  # # Save raw table
  write.table( x = output_table, 
               file = file.path( my_rp, "raw.txt"  ),
               sep = "\t", quote = F, row.names = T, col.names = T )
 
  write.table( x = geneset, 
               file = file.path( my_rp, "annotated_genes.txt"  ),
               sep = "\t", quote = F, row.names = F, col.names = F )
 
  # # Save pretty version
  pretty_table = output_table
  pretty_table$Genes  = NULL
  n_colors = length(unique(pretty_table$database)); 
  color_idx = pretty_table$database %>% factor(levels = desired_db, ordered = T) %>% as.integer
  my_cols = scales::hue_pal()(n_colors)[ color_idx ]
  theme_color_db = gridExtra::ttheme_minimal( core=list( bg_params = list( fill = my_cols ) ) )
  ggsave( gridExtra::tableGrob( d = pretty_table, theme = theme_color_db ), 
          file = file.path( my_rp, "color=database.pdf" ), 
          width = 15, height = N_ANNOT_PER_DB*length(desired_db) / 3, limitsize = F )
  
  return(output_table)
}

## ------------------------------------------------------------------------

#' Helper function. Check if a list of markers is compatible with a given Seurat object 
#' so that genes and cluster assignments are present in both `marker_info` and
#' the Seurat object `dge`.
#'
#' If `desired_cluster_order` is given, `are_compatible` checks that it is 
#' free of duplicates and it is a superset of the identity values occurring in other inputs.
are_compatible = function( dge, marker_info, ident.use ){
  atat( all( c("gene", "cluster") %in% names( marker_info ) ) )
  factor_flag = FALSE
  for( i in seq_along( marker_info ) ) {
    factor_flag = factor_flag || is.factor( marker_info[[i]] )
    marker_info[[i]] = as.character( marker_info[[i]] )
  }
  if( factor_flag ){ warning( "Factor columns detected in marker_info." ) }
  
  dge_ident =  as.character( FetchData( dge, ident.use )[[1]] )
  gene_compat  = all( marker_info$gene    %in% rownames( dge@data ) )
  ident_compat = all( marker_info$cluster %in% dge_ident )
  if( !gene_compat ){
    warning("marker_info$cluster has ID's not available in Seurat object")
  }
  if( !ident_compat ){
    warning("marker_info$gene has genes not available in Seurat object")
  }
  return( gene_compat && ident_compat )
}

#' Check if a list of cell types is good to be used to order genes and cells in a heatmap.
#' @details This means:
#' - no duplicates
#' - up to order, should be equal to union of table cluster labels and dge cluster labels
#' - no missing values
fix_cluster_order = function(  dge, marker_info, ident.use, desired_cluster_order = NULL ){
  if( is.null( desired_cluster_order ) ){
    warning( "No cluster order specified. Ordering clusters stupidly." )
    desired_cluster_order = union( Seurat::FetchData( dge, ident.use )[[1]], 
                                   marker_info$cluster )
  }
  atae( typeof( desired_cluster_order ), "character" )
  all_celltypes = union( Seurat::FetchData(dge, ident.use)[[1]], marker_info$cluster )
  missing = setdiff( all_celltypes, desired_cluster_order)
  extra   = setdiff( desired_cluster_order, all_celltypes)
  if( length( missing ) > 0 ){
    warning("Adding missing cell types to `desired_cluster_order` from Seurat object or marker table.")
    desired_cluster_order = c(desired_cluster_order, missing)
  }
  if( length( extra ) > 0 ){
    warning("Removing extraneous cell types from `desired_cluster_order`.")
    desired_cluster_order = intersect( desired_cluster_order, all_celltypes )
  }
  if( anyDuplicated( desired_cluster_order ) ){
    warning("Omitting duplicates in `desired_cluster_order`.")
    desired_cluster_order = unique( desired_cluster_order )
  }
  atat(  !any( is.na( desired_cluster_order ) ) )
  return( desired_cluster_order )
}


## ------------------------------------------------------------------------
#' Make a heatmap with one column for each cluster in `unique( Seurat::FetchData(dge, ident.use)[[1]])` and 
#' one row for every gene in `genes_in_order`. 
#' 
#' @param dge Seurat object
#' @param desired_cluster_order Levels of FetchData(dge, ident.use)[[1]], ordered how you want them to appear.
#' @param ident.use Variable to aggregate by.
#' @param labels "regular" (all labels), "stagger" (all labels, alternating left-right to fit more genes), or "none". 
#' @param aggregator Function to aggregate expression within a group of cells. Try mean or prop_nz.
#' @param normalize "row", "column", "both", or "none". Normalization function gets applied across the axis this specifies.
#' @param norm_fun Function to use for normalization. Try div_by_max or standardize.
#' @param genes_to_label A (small) subset of genes to label. If set, nullifies labels arg. 
#' @param main Title of plot.
#' @param return_type "plot" or "table". If "table", then instead of returning a heatmap, this 
#' returns the underlying matrix of normalized, cluster-aggregated values. If anything else, returns a ggplot.
#' @param ... Additional parameters passed to FetchData (such as use.raw = T).
#'
#' If the cluster's expression values are stored in `x`, then `aggregator(x)` gets (normalized and) plotted.
#' Optional parameter `desired_cluster_order` gets coerced to character. It should be a subset of 
#' `unique(Seurat::FetchData(dge, ident.use))` (no repeats).
#'
#' @export
#'
make_heatmap_for_table = function( dge, genes_in_order, 
                                   desired_cluster_order = NULL, 
                                   ident.use = "ident",
                                   labels = NULL, 
                                   aggregator = mean, 
                                   normalize = "row", 
                                   norm_fun = div_by_max,
                                   genes_to_label = NULL,
                                   main = "Genes aggregated by cluster",
                                   return_type = "plot", ... ){
  assertthat::assert_that(is.function(aggregator))
  assertthat::assert_that(is.function(norm_fun))
  if( !is.null( labels ) && !is.null(genes_to_label) ){
    warning("labels is ignored when genes_to_label is (are?) specified.\n")
    labels = "none"
  }
  
  # # Set up simple ident variable
  if(is.null(desired_cluster_order)){
    warning("No cell-type ordering given. Using arbitrary ordering.")
    desired_cluster_order = list(FetchData(dge, ident.use)[1, 1])
  }
  marker_info = data.frame( gene = genes_in_order, cluster = desired_cluster_order[[1]] )
  desired_cluster_order = fix_cluster_order( dge, marker_info, ident.use, desired_cluster_order )
  ident = FetchData(dge, ident.use) %>% 
    vectorize_preserving_rownames %>% 
    factor(levels = desired_cluster_order, ordered = T)
  ident = sort(ident)
  
  # # Sanitize input -- characters for genes, and no duplicate genes.
  genes_in_order = as.character( genes_in_order )
  if( anyDuplicated( genes_in_order ) ){
    warning( "Sorry, can't handle duplicate genes. Removing them." )
    genes_in_order = genes_in_order[ !duplicated( genes_in_order )]
  }
  
  # # Get cluster mean expression for each gene and row normalize
  logscale_expression = FetchDataZeroPad(dge, vars.all = genes_in_order, ...)[names( ident ), ]
  expression_by_cluster = aggregate_nice( x = logscale_expression, by = list( ident ), FUN = aggregator )
  # The net result of this conditional should be to preserve the dimensions.
  dim_old = dim(expression_by_cluster)
  if( normalize == "row" ){
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 2 )  
  } else if( normalize == "column" ){
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 1 ) %>% t 
  } else if( normalize == "both" ){
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 1 ) %>% t  
    expression_by_cluster = apply(X = expression_by_cluster, FUN = norm_fun, MARGIN = 2 ) 
  } else if( normalize != "none"){
    warning('normalize should be one of "row", "column", "both", or "none". Performing no normalization.')
    normalize = "none"
  } 
  atat(all(dim(expression_by_cluster) == dim_old))
  expression_by_cluster = t(expression_by_cluster) 
  
  
  # # Stop and return data if desired
  if( return_type == "table" ){ return(expression_by_cluster) }
  
  # # Form matrix in shape of heatmap and then melt into ggplot
  plot_df_wide = cbind( as.data.frame( expression_by_cluster ) , gene = rownames(expression_by_cluster))
  plot_df_wide$y = 1:nrow(plot_df_wide)
  plot_df_long = reshape2::melt( plot_df_wide, 
                                 id.vars = c("gene", "y"), 
                                 value.name = "RelLogExpr")
  plot_df_long$Cluster = factor( ordered = T, 
                                 as.character( plot_df_long$variable ),
                                 levels = desired_cluster_order )
  plot_df_long$gene = factor( ordered = T, 
                              as.character( plot_df_long$gene ),
                              levels = genes_in_order )

  p = ggplot( plot_df_long ) + ggtitle( main ) +
    geom_tile( aes(x = Cluster, y = gene, fill = RelLogExpr ) )
  p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Set default labeling mode according to number of genes
  if( is.null( labels ) ){ 
    if( length( genes_in_order ) < 30 ){
      labels = "regular"
    } else if ( length( genes_in_order ) < 80 ){
      labels = "stagger"
    } else {
      labels = "none"
    }
  }
  
  # # Remove labels if desired
  if( labels=="none"){
    p = p + theme(axis.ticks.y=element_blank(), axis.text.y = element_blank()) 
  }
  
  # # Make a DF with labeling info
  label_df = plot_df_long
  do_stagger = (labels == "stagger")
  label_df$horiz_pos = rep( c(-0.75*do_stagger, 0), length.out = nrow( plot_df_long ) )
  label_df %<>% extract(c("horiz_pos", "y", "gene")) 
  label_df = label_df[!duplicated(label_df$gene), ]
  invisible_label_row = label_df[1, ]
  invisible_label_row$gene = ""
  invisible_label_row$horiz_pos = -0.75 + min(label_df$horiz_pos)
  invisible_label_row$y = 1
  label_df = rbind(invisible_label_row, label_df)
  
  # # Add staggered labels with an invisible one farther out to make room (if desired)
  if( do_stagger ){
    p = p + theme(axis.ticks.y=element_blank(), axis.text.y = element_blank())
    p = p + geom_text(data = label_df, aes(x = horiz_pos, y = y, label = gene ))
  }
  
  # # Add only labels for genes_to_label (if desired)
  if( !is.null(genes_to_label) ){
    label_df$horiz_pos = 0.5
    label_df$horiz_pos[[1]] = -0.5
    label_df$gene[!(label_df$gene %in% genes_to_label)] = ""
    label_df = label_df[!duplicated(label_df$gene), ]
    label_layer = geom_text(data = label_df, aes(x = horiz_pos, y = y, label = gene )) 
      
    # # Attempt with ggrepel. Unsolved issue: labels migrate on top of heatmap.
    # label_layer = tryCatch( expr = { ggrepel::geom_label_repel(data = label_df, aes(x = horiz_pos, y = y, label = gene )) },
    #                         error = function(e){     geom_text(data = label_df, aes(x = horiz_pos, y = y, label = gene )) } )
    p = p + label_layer 

  }
  return( p )
}

#' Given a column ordering, produce a pleasing row ordering.
#'
#' @param X Set of data to order. Coerced to dataframe.
#' @param wrap Logical. Treat column n as a neighbor of column 1?
#' @param outlier_cutoff Numeric vector of length 1. See details.
#' @param REORDERFUN_outlier Function accepting a df X and returning a vector of indices. 
#' See details.
#' @param ... Extra args passed to REORDERFUN_outlier. 
#'
#' @return Returns a permutation of rownames(X).
#'
#' This function attempts to order rows for a heatmap to match a given column ordering. 
#' It attempts to create diagonal structure and it pays attention to which columns are adjacent.
#' The algorithm first orders by where the max value occurs: if row Alice peaks in column 1 and 
#' row Bob peaks in column 2, then Alice will precede Bob in the final ordering.
#' If Bob and Cassandra both peak in column 2, then the tiebreaker is column 3 minus column 1.
#'
#' @export
#'
OrderRowsForHeatmap = function( X, 
                                outlier_cutoff = 0.2, 
                                wrap = T,
                                REORDERFUN_outlier = function(x, ...) {
                                  x %>% dist %>% hclust %>% as.dendrogram %>% (stats::reorder)
                                }, ... ){
  X %<>% as.data.frame
  
  # Calculate max for each row.
  idx_max = apply(X, 1, which.max)
  idx_max_ctr   = cbind( 1:nrow(X), idx_max )
  idx_max_left  = cbind( 1:nrow(X), idx_max - 1 )
  idx_max_right = cbind( 1:nrow(X), idx_max + 1 )
  
  if( wrap ){
    # Alter indices: 0 becomes n, n+1 become
    has_neighbor_left  = rep( T, length( idx_max ) )
    has_neighbor_right = rep( T, length( idx_max ) )
    replace_a_with_b = function( x, a, b ){x[x==a]=b; return(x)}
    ncx = ncol(X)
    idx_max_left[, 2] %<>% replace_a_with_b( a = 0, b = ncx) %>% replace_a_with_b( a = ncx+1, b = 1)
    idx_max_right[, 2] %<>% replace_a_with_b( a = 0, b = ncx) %>% replace_a_with_b( a = ncx+1, b = 1)
  } else {
    # Remove indices referring to nonexistent entries
    has_neighbor_left  = ( idx_max!=1 )
    has_neighbor_right = ( idx_max!=ncol( X ) )
    idx_max_left  = idx_max_left[  has_neighbor_left, ]
    idx_max_right = idx_max_right[ has_neighbor_right, ]
  }
  
  # Calculate scores to order by
  major_score = idx_max
  minor_score = rep(0, length(major_score))
  minor_score[has_neighbor_left ] %<>% subtract( X[idx_max_left ] )
  minor_score[has_neighbor_right] %<>% add(      X[idx_max_right] )
  
  return(rownames(X)[order(major_score, minor_score, decreasing = F)])
  
  # Treatment of outlier rows is not yet implemented.
  # get_outlier_score = function( x, max_idx ){
  #   away_from_max = subset(seq_along(x), abs( seq_along(x) - max_idx ) > 1 ) 
  #   sum(x[away_from_max]) / sum(x)
  # }
  
}


## ------------------------------------------------------------------------
#' @export
screen_receptor_ligand = function( is_expressed, results_path ){
  
  # # Get receptor-ligand pairs; annotate with tissues expressed; save
  ramilowski = get_ramilowski()
  ramilowski$Ligand.ApprovedSymbol = NULL
  ramilowski$Receptor.ApprovedSymbol = NULL
  ramilowski = subset( ramilowski, ligand_mouse %in% rownames(is_expressed) )
  ramilowski = subset( ramilowski, receptor_mouse %in% rownames(is_expressed) )
  ramilowski$ligand_cell_types = 
    is_expressed[ramilowski$ligand_mouse, ] %>% 
    apply( 1, which ) %>% 
    sapply( names ) %>% 
    sapply( paste, collapse = "_")
  ramilowski$receptor_cell_types = 
    is_expressed[ramilowski$receptor_mouse, ] %>% 
    apply( 1, which ) %>% 
    sapply( names ) %>% 
    sapply( paste, collapse = "_")
  write.table( ramilowski, file =  file.path( results_path, "Receptor_ligand_all.txt" ), 
               sep = "\t", row.names = F, col.names = T, quote = F )

  absent = union( ramilowski$ligand_mouse, ramilowski$receptor_mouse ) %>% setdiff( rownames( is_expressed ) )
  if( length( absent ) > 0 ){
    zeropad = matrix(F, ncol = is_expressed, nrow = length( absent ), 
                     dimnames = list( gene = absent, 
                                      celltype = colnames(is_expressed)) )
    is_expressed %<>% rbind( zeropad )
  }
  
  # # Get lists of receptors and ligands for each tissue pairing
  num_unique_ligands = matrix( 0, nrow = 3, ncol = 3, 
                               dimnames = list( lig_expr_tissue = c("BLD", "MES", "TEC"),
                                                rec_expr_tissue = c("BLD", "MES", "TEC") ) )
  dir.create.nice( file.path( results_path, "ligand_lists" ) )
  dir.create.nice( file.path( results_path, "receptor_lists" ) )
  for( lig_tissue in rownames(num_unique_ligands)){
    for( rec_tissue in colnames(num_unique_ligands)){
      eligible_subset = subset( ramilowski, 
                                is_expressed[ligand_mouse,   lig_tissue] & 
                                  is_expressed[receptor_mouse, rec_tissue] )
      num_unique_ligands[lig_tissue, rec_tissue] = length( unique( eligible_subset$ligand_mouse ) )
      write.table( unique(eligible_subset$ligand_mouse  ), row.names = F, col.names = F, quote = F,
                   paste0( results_path, "/ligand_lists/"  , lig_tissue, "_to_", rec_tissue, ".txt") )
      write.table( unique(eligible_subset$receptor_mouse), row.names = F, col.names = F, quote = F,
                   paste0( results_path, "/receptor_lists/", lig_tissue, "_to_", rec_tissue, ".txt") )
      
    }  
  }
  
  # Background lists composed of everything that got successfully converted to mouse ortholog
  # For pathway analysis with a background list
  ramilowski_orig = read.table( file.path( PATH_TO_TABLES, "LigandReceptor_Ramilowski2015_mouse.txt" ), 
                                header = T, sep="\t", stringsAsFactors = F )
  write.table( ramilowski_orig$ligand_mouse   %>% unique, 
               paste0( results_path, "/ligand_lists/background.txt"),
               row.names = F, col.names = F, quote = F)
  write.table( ramilowski_orig$receptor_mouse %>% unique, 
               paste0( results_path, "/receptor_lists/background.txt"),
               row.names = F, col.names = F, quote = F)
  
  
  # # Save summary to file
  # # Sink helps get the full dimnames
  sink( file.path( results_path, "num_unique_ligands.txt" ) )
  {
    print( num_unique_ligands )
  }
  sink()
  return()
}

## ------------------------------------------------------------------------

#' Quickly explore many parameter settings
#'
#' @param dge Seurat object.
#' @param results_path 
#' @param all_params Dataframe specifying multiple runs of the pipeline. Must have these columns: 
#' \itemize{
#' \item "latent_dimension": integer. Latent dimension to use in PCA or CCA.
#' \item "clust_method": Character. Each entry should be a valid input for the 'method' arg of 'cluster_wrapper'. (Currently 'SNN' or 'DBSCAN').
#' \item "clust_granularities_as_string": Character. Comma-separated list of numbers to be used for G.use (if DBSCAN) or res (if SNN). Allows multiple clusterings based on the same PC list. 
#' \item "excess_var_cutoff": Numeric. Params for variable gene selection. Units are in multiples of the local median coefficient of variation, i.e. 0.8 means a gene is included if it exceeds 80% of the local median of the CV, binning based on expression level.
#' \item "log_expr_cutoff": Numeric. Genes below this expression level are excluded. Done by averaging the '@data' slot (so, log1p-scale expression).
#' }
#' @param whitelist Everything not on this list is excluded.
#' @param blacklist Genes to exclude. Typically nothing (default), but could be cell cycle genes (try 
#' @param regress_out Covariates to regress out via Seurat::RegressOut. . 
#' \code{ thymusatlastools2::get_macosko_cc_genes() \%>\% Reduce(f=union) } or ribosomal 
#' subunit genes. Make sure you get the species right!
#' @param plot_each_iteration A list of featureplots to show at each iteration.
#' @param do_cca @param align_embeddings If do_cca, do CCA instead of PCA. If align_embeddings, align results using DTW.
#' @param group.by @param group1 @param group2 @param ... Passed to RunCCA. 
#'
#' @export

explore_embeddings = function( dge, results_path, all_params, 
                               whitelist = AvailableData(dge), 
                               blacklist = c(), 
                               regress_out = c(),
                               plot_each_iteration = "orig.ident",  
                               do_cca = F, align_embeddings = F, 
                               group.by = NULL, group1 = NULL, group2 = NULL,
                               ... ){
  # Check to see all required inputs are present
  atat(is.data.frame( all_params ) )
  if (is.null(all_params$latent_dimension)){
    all_params$latent_dimension = all_params$num_pc
  }
  required_params = c( "latent_dimension", "clust_granularities_as_string",
                       "excess_var_cutoff","log_expr_cutoff" )
  atat( all( required_params %in% names( all_params ) ) )
  
  if(do_cca){
    assertthat::assert_that( !is.null(group.by) & !is.null(group1) & !is.null(group2) ) 
  }
  
  # # Record the things you're gonna try
  dir.create.nice( results_path )
  write.table( all_params, file = file.path(results_path, "params_to_try.txt"), 
               quote = F, sep = "\t", row.names = F)
  
  # # Try all the things
  prev_param_row = NULL
  for(i in rownames( all_params ) ){
    param_row = all_params[i, ]; names(param_row) = names(all_params)
    print("Trying these settings:")
    print(param_row)
    rp_mini = file.path(results_path, collapse_by_name( all_params[i,] ))
    dir.create.nice(rp_mini)
    
    # # select genes, confining to whitelist and excluding blacklist.
    dge = var_gene_select( dge, results_path = rp_mini,
                           excess_var_cutoff   = param_row[["excess_var_cutoff"]],
                           log_expr_cutoff     = param_row[["log_expr_cutoff"]] )
    dge@var.genes %<>% intersect( whitelist )
    dge@var.genes %<>% setdiff( blacklist )
    
    # # Scale genes
    dge %<>% FillNA
    if( length( regress_out ) > 0 ){
      dge %<>% Seurat::ScaleData(vars.to.regress = regress_out, check.for.norm = F, genes.use = dge@var.genes)
    } else {
      dge %<>% Seurat::ScaleData( check.for.norm = F, genes.use = dge@var.genes  )
    }
    
    # Reduce dimension and cluster cells
    if( do_cca ){
      dge = Seurat::RunCCA(dge, pcs.compute = param_row[["latent_dimension"]], do.print = F, 
                           group.by = group.by, group1 = group1, group2 = group2, ... ) 
      reduction.use = "cca"
    } else {
      dge = Seurat::RunPCA(dge, pcs.compute = param_row[["latent_dimension"]], do.print = F  ) 
      reduction.use = "pca"
    }
  
    if( align_embeddings ){
      dge = Seurat::AlignSubspace(dge, 
                                  reduction.type = reduction.use, 
                                  show.plots = F,
                                  dims.align = 1:param_row[["latent_dimension"]], 
                                  grouping.var = group.by ) 
      reduction.use = paste0(reduction.use, ".aligned")
    }
    
    dge %<>% Seurat::RunTSNE(reduction.use = reduction.use, 
                             dims.use = 1:param_row[["latent_dimension"]], 
                             check_duplicates = FALSE) 
    dge %<>% Seurat::RunUMAP(reduction.use = reduction.use, 
                             dims.use = 1:param_row[["latent_dimension"]] ) 
    G = param_row[["clust_granularities_as_string"]] %>% strsplit(",") %>% extract2(1) %>% as.numeric
    dge %<>% FindClusters( reduction.type = reduction.use,
                           dims.use = 1:param_row[["latent_dimension"]],
                           resolution = G, force.recalc = T )
    
    # # save plots and summaries 
    saveRDS( dge, file.path( rp_mini, "dge.data") ) 
    if(length(plot_each_iteration ) > 0 ){
      save_feature_plots(dge, results_path = rp_mini,
                         gene_list = plot_each_iteration, gene_list_name = "summary_plots",
                         types = "PDF")
      save_feature_plots(dge, results_path = rp_mini, axes = c("UMAP1", "UMAP2"),
                         gene_list = plot_each_iteration, gene_list_name = "summary_plots",
                         types = "PDF")
    }
    prev_param_row = param_row
  }
  return(dge)
}


## ------------------------------------------------------------------------
#' Create an information-dense metadata variable: clusters annotated with marker genes.
#'
#' @param dge Seurat object
#' @param results_path Where to put the resulting plot. If NULL, no plot is saved.
#' @param genes.use 
#' @param gene_list_name Goes into name of resulting plot.
#' @param n Number of genes to label each cluster with
#' @param ident.use Clusters to annotate
#' @param ... Additional parameters for annot_ident_plot
#'
#' @export
#'
InstantOverview = function( dge, results_path = NULL,
                            genes.use, gene_list_name = deparse(substitute(genes.use)),
                            n = 5, ident.use = "ident", include_size = T, ... ){
  # Parse input
  ident = FetchData(dge, ident.use)[[1]]
  if(!is.factor(ident)){
    warning("To guarantee color matched plots, use a factor variable for ident.use.")
    ident_levels = ident %>% unique %>% sort
  } else {
    ident_levels = levels( ident )
  }
  
  # Create size string (empty if not desired)
  if(include_size){
    cluster_sizes = table(ident)[ident_levels] %>% paste0( " (", ., ") ")
  } else {
    cluster_sizes = rep("", length( ident_levels ) )
  }
  
  # Create annotated labels
  dge %<>% AddClusterIndicators(ident.use = ident.use)
  ident_plus_gene_levels = ident_levels; names(ident_plus_gene_levels) = ident_levels
  for( i in seq_along( ident_levels )){
    indicator_name = paste0( ident.use, ident_levels[[i]])
    label_genes = dge %>% get_similar_genes( markers = indicator_name, n = n, genes.use = genes.use )
    ident_plus_gene_levels[[i]] %<>% paste0( cluster_sizes[[i]], "\n", paste0( label_genes, collapse = "\n") ) 
  }
  
  # Add annotated levels to Seurat object
  translate_levels = function(x) ident_plus_gene_levels[as.character( x )]
  dge@meta.data$ident_plus_markers = translate_levels( ident ) %>% 
    factor( ordered = T, levels = translate_levels( ident_levels ) )
  
  # Make plots
  annot_ident_plot( dge, 
                    results_path,
                    figname = "annot_ident_plot",
                    ident.use = ident.use, 
                    ... )
  annot_ident_plot( dge,
                    results_path,
                    figname = "annot_ident_plot" %>% paste0("_", gene_list_name), 
                    ident.use = "ident_plus_markers", 
                    ... )
  return( dge )
  
}

