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
  plot_df_long$Cluster = factor( as.character( plot_df_long$variable ),
                                 levels = desired_cluster_order )
  plot_df_long$gene = factor( as.character( plot_df_long$gene ),
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


#' Return a cell ordering where predefined clusters remain contiguous.
#' 
OrderCellsWithinClusters = function( dge, ident.use, coords.use = paste0("PC", 1:30), cluster_order = NULL ){
  
  # The implementation scheme here is to apply a series of permutations,
  # first to order the cells by cluster and then to order the cells within each cluster.
  # At the end, I can produce a single reordering equivalent to the composition of the
  # steps described above. To do that, I'll keep track of the original position.
  
  # Set up df with ident and original position
  clusters_and_orig_pos = Seurat::FetchData(dge, ident.use)
  clusters_and_orig_pos[["orig_pos"]] = seq_along(clusters_and_orig_pos[[ident.use]])
  clusters_and_orig_pos[[ident.use]] %<>% as.character
    
  # Order cells by cluster
  if(!is.null(cluster_order)){
    expected_levels = unique(clusters_and_orig_pos[[ident.use]])
    assertthat::are_equal(sort(cluster_order), sort(expected_levels))
    clusters_and_orig_pos[[ident.use]] %<>% factor( ordered = T, levels = cluster_order )
  }
  clusters_and_orig_pos = clusters_and_orig_pos[order(clusters_and_orig_pos[[ident.use]]), ]
 
  # Reorder within each cluster
  for( cluster in unique( clusters_and_orig_pos[[ident.use]] ) ){
    this_cluster_idx_current = which( clusters_and_orig_pos[[ident.use]] == cluster )
    this_cluster_idx_orig = clusters_and_orig_pos[this_cluster_idx_current, "orig_pos"]
    this_cluster_coords = Seurat::FetchData(dge, coords.use)[this_cluster_idx_orig, ]
    this_cluster_reordering = this_cluster_coords %>% dist %>% (fastcluster::hclust) %>% extract2("order")
    this_cluster_idx_new = this_cluster_idx_current[this_cluster_reordering]
    clusters_and_orig_pos[this_cluster_idx_current, ] = clusters_and_orig_pos[this_cluster_idx_new, ] 
  }
  
  # Return a reordering that works in one shot.
  return( clusters_and_orig_pos[["orig_pos"]] )
}

#' Save a big PDF file to `<results_path>/<main>.pdf` containing a heatmap of gene expression levels.
#' 
#' @param dge a Seurat object
#' @param results_path should be a character such that `dir.exists( results_path )`
#' @param genes.use Genes to include in the heatmap. Character vector.
#' @param genes.preview If this has length > 0, function will produce a quick preview containing these genes.
#' @param norm_fun Applied to each gene. 
#' @param num_pc Number of PC's to use when reordering cells
#' @param main Figure title and name of saved file.
#' @param ident.use Used to set up the ColSideColors in heatmap.2.
#' @param width @param height Dimensions of saved plot in inches.
#' @param ... Passed on to heatmap.2.
#' @param cluster_order What order do you want the clusters in? 
#' @param cluster_colors List of colors. Should be named with levels of the FetchData(dge, ident.use).  
#' Used to set up the ColSideColors in heatmap.2. If you don't want any, put NULL. 
#' If NA, defaults to a sensible value.
#' @param col a vector of color names for use in the main body of the heatmap.
#' @details Each column is a cell and each row is a gene. Each gene is rescaled so that its peak expression is 1.
#' This facilitates comparison within genes and across cells, though it's bad for comparison across genes.
#'
#'@export
#'
DoHeatmapFast = function( dge, results_path, 
                          genes.use = dge@var.genes, 
                          genes.preview = c(),
                          genes.label = c(),
                          main = "heatmap_cellwise",
                          ident.use = "ident",
                          norm_fun = div_by_max, 
                          cluster_colors = NA,
                          cluster_order = NULL,
                          num_pc = 20,
                          col = colorRampPalette(c("khaki1", "red"))(20), 
                          width = 10, height = 10,
                          ... ){
  if( length(genes.preview) == 1 | length(genes.use) == 1 ){
    stop("Sorry, this function is not guaranteed to work with only one gene for genes.use or genes.preview.")
  }
  if( is.na(cluster_colors)){
    cluster_names = unique(Seurat::FetchData(dge, ident.use)[[1]])
    cluster_colors = scales::hue_pal()(length(cluster_names))
    names(cluster_colors) = cluster_names
  }

  cat("Projecting PCA...\n")
  dge %<>% (Seurat::ProjectPCA)( do.print = F, pcs.store = num_pc )
  #save memory; these aren't used downstream
  dge@raw.data = matrix()
  dge@scale.data = matrix()
  
  # cat("Ordering genes...\n")
  # genes.use = dge@pca.x.full[genes.use, ] %>% dist %>% (fastcluster::hclust) %>% extract2("order")
  # genes.use = dge@var.genes[genes.use]
  
  # reorder cells and set sparse row labels
  cat("Ordering cells...\n")
  cell_order = OrderCellsWithinClusters( dge, ident.use = ident.use, coords.use = paste0("PC", 1:num_pc), 
                                         cluster_order = cluster_order )
  
  # sweep out max expression level and set colorscale
  cat("Normalizing expression...\n")
  norm_expr = t( apply(X = dge@data[genes.use, cell_order], FUN = norm_fun, MARGIN = 1) )
  gene_labels = genes.use
  names(gene_labels) = genes.use
  gene_labels[ !(gene_labels %in% genes.label) ] = ""

  if( length(genes.preview)>0 ){
    cat("Making preview... \n")
    fname = paste0( main, "_PREVIEW.pdf" )
    pdf( file.path( results_path, fname ), width = width, height = height )
    genes.preview = rownames(norm_expr) [rownames(norm_expr) %in% genes.preview] # like intersect(), but preserves order
    preview_cells = sample( 1:ncol( norm_expr ), min( 300, ncol( norm_expr ) ) ) %>% sort
    colorbar_subset = cluster_colors[ as.character(Seurat::FetchData(dge, ident.use)[cell_order[preview_cells], 1]) ]
    gplots::heatmap.2( norm_expr[genes.preview, preview_cells], 
                       Rowv = F, 
                       Colv = F, 
                       dendrogram = "none",
                       symm = F, 
                       scale = "none", 
                       col = col,
                       trace = "none",
                       xlab = "Cells", labCol = "",
                       ylab = "Genes", labRow = gene_labels[genes.preview],
                       ColSideColors = colorbar_subset, 
                       ... )
    dev.off()
  }

  
  cat("Making full heatmap...\n")
  fname = paste0( main, ".pdf" )
  pdf( file.path( results_path, fname ), width = width, height = height  )
  gplots::heatmap.2( norm_expr, 
                     Rowv = F, 
                     Colv = F, 
                     dendrogram = "none",
                     symm = F, 
                     scale = "none", 
                     col = col,
                     trace = "none",
                     xlab = "Cells", labCol = "",
                     ylab = "Genes", labRow = gene_labels,
                     ColSideColors = cluster_colors[ as.character(Seurat::FetchData(dge, ident.use)[cell_order, 1]) ],
                     ... )
  dev.off()
  cat("Done.\n")
}


## ------------------------------------------------------------------------
#' Add a colourbar to a ggplot object, fixing the fill to allow a separate colorscale. 
#' 
gg_add_colorbar = function( plot, 
                            x, 
                            my_labels,
                            col = NULL,
                            width_ = 0.1, 
                            position = -0.1, 
                            is_horiz = T ){
  my_labels %<>% make.names
  unique_labels = my_labels %>% unique
  if( is.null( col ) ){
    col = setNames( 
      scales::hue_pal()( my_labels %>% unique %>% length ),
      my_labels %>% unique %>% make.names 
    )
  } else {
    if( !all( unique_labels %in% names( col ) ) ){
      warning( "col should be named with some permutation of unique(make.names(my_labels)).\n" )
    }
  }
  
  colourbar_obs = data.frame( my_labels = my_labels, x = x, position = position - 0.1*replace_with_int_rank(my_labels) )
  n = nrow(colourbar_obs)
  for( cluster in unique_labels ){
    if( is_horiz ){
      plot = plot + 
        geom_tile( data = subset( colourbar_obs, my_labels == cluster ),
                   mapping = aes( x = x, y = position ), 
                   height = width_,
                   fill = col[cluster] ) 
      # geom_density( data = subset( colourbar_obs, my_labels == cluster ),
      #               mapping = aes( x = x, y = -n*0.2*..density.. ),
      #               fill = col[cluster], alpha = 0.5 )
    } else {
      plot = plot +
        geom_tile( data = subset( colourbar_obs, my_labels == cluster ),
                   mapping = aes( y = x, x = position),
                   width = width_,
                   fill = col[cluster] ) 
    }
  }       
  return(plot)
}

