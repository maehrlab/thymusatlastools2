## ------------------------------------------------------------------------


#' Get the order of the tips from an object of class phylo from the APE package.
#'
extract_tip_order = function( ape_phylo_object ){
  is_tip = ape_phylo_object$edge[,2] <= length(ape_phylo_object$tip.label)
  ordered_tips = ape_phylo_object$edge[is_tip, 2]
  ape_phylo_object$tip.label[ordered_tips]
}

#' Re-number clusters to be compatible with dendrogram.
#'
#' @export
#' 
RenameByDendro = function ( object ){
  if( is.null( object@cluster.tree ) ){
    object %<>% BuildClusterTree
  }
  object %<>% BuildClusterTree
  auto_names = extract_tip_order(object@cluster.tree[[1]]) %>% as.character
  phylo_names = paste0(1:length(auto_names))
  converter = setNames(phylo_names, auto_names)
  object %<>% SetIdent( ident.use = converter[as.character(object@ident)] )
  object@cluster.tree[[1]]$tip.label %<>% extract(converter, .)
  return( object )
}


#' Remove missing values from the metadata, issuing a warning if changes are made. 
#'
#' @param object Seurat object
#'
#' @export
#'
FillNA = function( object, filler = "NA" ){
  varnames = object@meta.data %>% names
  missing_list = c()
  na2filler = function(x){
    x[is.na(x)] = filler
    return(x)
  }
  for( var in varnames ){
    if( any( is.na( object@meta.data[[var]] ) ) ){
      missing_list %<>% c(var)
      object@meta.data[[var]] %<>% na2filler
    }
  }
  if( length( missing_list ) > 0 ){
    warning( paste0( "Missing values found in these variables: \n", 
                     paste0( missing_list, collapse = "\n" ),
                     "\nReplacing with ", filler, " .\n\n" ) )
  }
  return( object )
}


#' Get all available variable names (genes, identity classes, PCA embeddings, etc)
#'
#' @param object Seurat object
#' @return Returns a character vector of all eligible inputs to the `vars.all` argument of `FetchData`.
#' @export
#'
AvailableData = function( object ){
  get_dr_names = function(x) colnames(x@cell.embeddings)
  dr_names = object@dr %>% lapply( get_dr_names ) %>% Reduce(f=union)
  available_categorized = list( metadata = names( object@meta.data ),
                                dr_names,
                                genes = rownames( object@data ),
                                ident = "ident" )
  assay_data = lapply( object@assay, function(x) rownames(x@data) )
  available_categorized %<>% c(assay_data)
  return( Reduce( f = union, x = available_categorized ) )
}

#' Sanitize gene names via `make.names`
#'
#' @export 
#'
SanitizeGenes = function( dge ){
  mn_null_tolerant = function(x) {
    if(is.null(x)){
      return(NULL)
    } else {
      return(make.names(x))
    }
  }
  rownames( dge@raw.data )   %<>% mn_null_tolerant
  rownames( dge@data )       %<>% mn_null_tolerant
  rownames( dge@scale.data ) %<>% mn_null_tolerant
  names(    dge@var.genes )  %<>% mn_null_tolerant
  return( dge )
}


#' FetchData but with zeroes for unavailable genes
#'
#' @export
#' @param dge Seurat object
#' @param vars.all List of all variables to fetch. Missing entries are ignored.
#' @param ... Other arguments to pass to FetchData
#'
#' @details This function is stupid: if you ask for "PC1" and it's not available,
#' it will think you're asking for a non-expressed gene, so it will return zeroes.
FetchDataZeroPad = function( dge, vars.all, warn = T, ... ){
  vars.all = vars.all[complete.cases(vars.all)]
  avail = intersect( vars.all, AvailableData( dge ) )
  unavail = setdiff( vars.all, AvailableData( dge ) )
  if(warn && length(unavail) > 0){
    warning("Some variables are not available. Returning zeroes.\n")
  }
  to_return  = FetchData( dge,  avail, ... ) 
  pad = as.data.frame( matrix(0,           
                              nrow = nrow( to_return ), 
                              ncol = length( unavail ),
                              dimnames = list( rownames( to_return ),               
                                               unavail) ) )
  to_return = cbind( to_return, pad )
  assertthat::are_equal( sort( vars.all ),   sort( colnames( to_return ) ) )   
  to_return = to_return[, vars.all, drop = F]
  assertthat::assert_that( is.data.frame( to_return ) )
  return( to_return )
}


#' Subset data flexibly from a Seurat object.
#'
#' @param dge Seurat object
#' @param vars.use Variables to fetch for use in `predicate`.
#' @param predicate String to be parsed into an R expression and evaluated as an arg to base::subset.
#' @param preserve_raw If TRUE, then it will not exclude any cells from the @raw.data slot.
#' By default, it leaves cells in @raw.data only if they satisfy the given predicate.
#' @param print_numbers If true, save a text file with cellcounts (initial, removed, retained).
#' @param show_on_tsne If true, save a tSNE showing which cells are kept.
#' @param results_path Where to save the tSNE.
#' @param ... Extra params passed to tsne_colored.
#' @details Calls FetchData, subsets it as a dataframe using base::subset, and 
#' subsets the Seurat object correspondingly (using the df rownames.)
#'
#' @export
#'
SubsetDataFlex = function( dge, vars.use, predicate, preserve_raw = FALSE, 
                           results_path = NULL, 
                           print_numbers = !is.null(results_path),
                           show_on_tsne = !is.null(results_path), ... ){
  if( typeof(predicate) != "character"){
    print("predicate should be a character vector. It will be parsed into `subset` as an R expression.")
  }

  df = FetchData(dge, vars.use) %>% as.data.frame
  parent_env = parent.frame()
  cu = df %>% subset(eval( expr = parse(text=predicate), envir = df, enclos = parent_env )) %>% rownames
  if( show_on_tsne ){
    assertthat::assert_that(!is.null(results_path))
    dge@meta.data$included = (rownames(dge@meta.data) %in% cu) %>% as.character
    tsne_colored( dge, results_path, colour = "included", ... )
  }
  
  if( print_numbers ){
    inc = (rownames(dge@meta.data) %in% cu) 
    assertthat::assert_that(!is.null(results_path))
    X = data.frame(
      category = c("initial", 
                   "removed", 
                   "retained"),
      count = c(length(inc), 
                sum(!inc), 
                sum(inc))
    )
    write.csv(X, file.path(results_path, "cellcounts.csv"))
  }
  
  dge = SubsetData(dge, cells.use = cu)
  if( !preserve_raw ){
    dge@raw.data = deseuratify_raw_data(dge) 
  }
  return( dge )
}

#' Merge two Seurat objects.
#'
#' By default, preserves any column of @meta.data shared between both objects.  
#' You can also specify what variables to keep. They will be added to meta.data in 
#' the output, warning and padding with zeroes if either object lacks any var in vars.keep.
#'  
#' @export
#' 
SeuratMerge = function( dge1, dge2, preserve_ident = F, 
                        vars.keep = intersect(names(dge1@meta.data), 
                                              names(dge2@meta.data) ) ){
  # Save on memory
  dge1@scale.data = matrix()
  dge2@scale.data = matrix()
  gc()
  
  # Do merging
  dge_all = list( dge1 = deseuratify_raw_data( dge1 ), 
                  dge2 = deseuratify_raw_data( dge2 ) ) %>%
    dge_merge_list %>% seuratify_thy_data 
  characterize_factor = function(x){ if(is.factor(x)) as.character(x) else x }
  all_factors_to_char = function(X) data.frame(lapply(X, characterize_factor), stringsAsFactors=FALSE)
  
  if( length(vars.keep) > 0 ){
    preserved_metadata = rbind( FetchDataZeroPad( dge1, vars.keep ) %>% all_factors_to_char, 
                                FetchDataZeroPad( dge2, vars.keep ) %>% all_factors_to_char )
    preserved_metadata %<>% as.data.frame
    rownames(preserved_metadata) = c( rownames(dge1@meta.data),rownames(dge2@meta.data)) 
    dge_all %<>% AddMetaData(preserved_metadata)
  }
  
  if(preserve_ident){
    new_ident = c( characterize_factor(dge1@ident), 
                   characterize_factor(dge2@ident) )
    names(new_ident) = c(names(dge1@ident), names(dge2@ident))
    dge_all %<>% SetIdent(new_ident, cells.use = names(new_ident) )
  }
  
  return(dge_all)
}



## ------------------------------------------------------------------------
#' Create indicator vectors for a categorical variable.
#'
#' @param dge Seurat object
#' @param results_path If not null, save plots here. 
#' @param ident.use Categorical vector to create indicators for
#'
#' @export
#'
AddClusterIndicators = function( dge, results_path = NULL, ident.use = "ident", ... ){
  ident_df = data.frame(FetchData(dge, ident.use) )
  colnames( ident_df ) = ident.use
  my_form = paste0("~", ident.use, " + 0", collapse = "") %>% as.formula
  cluster_indicators = model.matrix( my_form, ident_df ) %>% as.data.frame
  dge %<>% AddMetaData(cluster_indicators)
  if(!is.null( results_path )){
    save_feature_plots( dge, results_path, 
                        gene_list = colnames(cluster_indicators), 
                        ...)
  }
  return( dge )
}

## ------------------------------------------------------------------------

#' Make a FACS-like plot from a single-cell rna-seq dataset.
#'
#' @param dge Seurat object
#' @param gene1 Horizontal axis on plot mimics this gene. Character, usually length 1 but possibly longer.
#' @param gene2 Vertical axis on plot mimics this gene. Character, usually length 1 but possibly longer. 
#' @param genesets_predetermined If FALSE, plot the sum of many genes similar to gene1 instead of gene1 alone (same 
#' for gene2). See ?get_similar_genes. If TRUE, plot the sum of only the genes given.
#' @param num_genes_add Each axis shows a simple sum of similar genes. This is how many (before removing overlap). Integer.
#' @param return_val If "all", returns a list with several internal calculations revealed.
#' If "plot", returns just a ggplot object. If "seurat", returns a Seurat object with gene scores added. 
#' @param cutoffs If given, divide plot into four quadrants and annotate with percentages. Numeric vector of length 2.
#' @param dge_reference Seurat object. This function relies on gene-gene correlation. If your dataset is perturbed in a way 
#' that would substantially alter gene-gene correlations, for example if different time points are present or certain 
#' cell types are mostly depleted, you can feed in a reference dge, and TACS will choose axes based on the reference data.

#' @param density If TRUE, plot contours instead of points.
#' @param ... Extra params for stat_density2d.
#'
#' This function is based on a simple scheme: choose genes similar to the ones specified 
#' and average them to reduce the noise. 
#'
#' @export 
#'
TACS = function( dge, gene1, gene2, cutoffs = NULL, 
                 return_val = "plot", density = F, 
                 facet_by = NULL, 
                 include_panel_with_all = FALSE, 
                 facet_levels = 
                   FetchData(dge, facet_by)[[1]] %>% factor %>% levels %>% 
                   c(rep("all", include_panel_with_all), .),
                 col = stats::setNames( scales::hue_pal()( length( facet_levels ) - include_panel_with_all ), 
                                        facet_levels[ ( include_panel_with_all + 1 ) : length( facet_levels )] ),
                 num_genes_add = 100, genesets_predetermined = F, dge_reference = dge, ... ){
  # Get gene sets to average
  if(genesets_predetermined){
    g1_similar = gene1
    g2_similar = gene2
  } else {
    g1_similar = get_similar_genes(dge_reference, gene1, num_genes_add) %>% c( gene1, . )
    g2_similar = get_similar_genes(dge_reference, gene2, num_genes_add) %>% c( gene2, . ) 
    shared = intersect(g1_similar, g2_similar)
    g1_similar %<>% setdiff(shared)
    g2_similar %<>% setdiff(shared)
  }
  
  # Average gene sets to get scores
  g1_score = rowMeans(FetchDataZeroPad(dge, g1_similar))
  g2_score = rowMeans(FetchDataZeroPad(dge, g2_similar))
  g1_score_name = paste0(gene1[1], "_score")
  g2_score_name = paste0(gene2[1], "_score")
  
  #Add scores as metadata. Extract with faceting var into plotting data.
  dge %<>% AddMetaData(g1_score, col.name = g1_score_name)
  dge %<>% AddMetaData(g2_score, col.name = g2_score_name)
  plot_df = FetchData(dge, c(g1_score_name, g2_score_name, facet_by))
  
  # Augment data to form extra panel with everything
  if( include_panel_with_all ){
    plot_df_all = plot_df
    plot_df_all[[facet_by]] = "all"
    plot_df = rbind(plot_df, plot_df_all)
    col = c(col, "all"="black")
  } 
  # Prepare to facet
  if(!is.null(facet_by)) {
    plot_df[[facet_by]] %<>% factor(levels = facet_levels, ordered = T) %>% droplevels
  }
  
  # Form plot
  p = ggplot(plot_df) 
  if(density){ 
    p = p + stat_density2d( aes_string( x = g1_score_name, y = g2_score_name, 
                                        colour = facet_by, alpha = "..level.." ), bins = 50 ) +
      scale_alpha_continuous( range = c(0.4, 1) ) + 
      scale_color_manual(values = col)
  } else {
    p = p + geom_point( aes_string( x=g1_score_name, y=g2_score_name ) ) 
  }
  p = p + expand_limits(y=0, x=0)
  # Facet if desired
  if(!is.null(facet_by)) {
    p = p + facet_wrap(as.formula(paste0("~", facet_by)))
  }
  
  # Add quadrants and percentages
  if( !is.null(cutoffs)){
    p %<>% add_quadrants(g1_score_name = g1_score_name,
                         g2_score_name = g2_score_name, 
                         cutoffs = cutoffs,
                         facet_by = facet_by)
  } 
  
  # Return everything or just a plot or just a seurat object
  if( return_val == "all" ){
    return( list( plot = p, 
                  dge = dge, 
                  genes = list( gene1 = gene1, gene2 = gene2 ),
                  score_names = c( g1_score_name, g2_score_name ), 
                  genesets = list( g1_similar, g2_similar ),
                  cutoffs = cutoffs,
                  plot_df = plot_df ) )
  } else if( return_val == "seurat" ){
    return(dge)
  } else if( return_val == "plot" ){
    return( p )
  } else {
    stop(" return_val should be 'all', 'seurat', or 'plot'. ")
  }
}

#' Split a scatterplot (or similar) into quadrants and label percentages in each quadrant. 
#'
#' @param dge Seurat object
#' @param cutoffs numeric vector of length 2. Where to delineate the quadrants.
#' @param facet_by optional facetting variable. Percents are calculated separately for each facet.
#'
#' This is a helper for TACS, but it's exported in case it becomes useful.
#'
#' @export 
#'
add_quadrants = function(p, g1_score_name, g2_score_name, cutoffs, facet_by = NULL){
  
  # Calculate percentages
  p = p + geom_vline(data = data.frame(xint=cutoffs[1]), aes(xintercept=xint))
  p = p + geom_hline(data = data.frame(yint=cutoffs[2]), aes(yintercept=yint))
  percentages = p$data[c(g1_score_name, g2_score_name, facet_by)]
  percentages[[g1_score_name]] %<>% is_greater_than(cutoffs[1])
  percentages[[g2_score_name]] %<>% is_greater_than(cutoffs[2])
  percentages %<>% table %>% (reshape2::melt)
  if(!is.null(facet_by)) {
    percentages = percentages[order(percentages[[facet_by]]), ]
    for( facet_level in unique(p$data[[facet_by]])){
      percentages[percentages[[facet_by]] == facet_level, "value"] %<>% percentify()
    } 
  } else {
    percentages$value %<>% percentify()
  }
  
  # Form annotation DF with correct facet and attempt sensible placement of percentages
  for( i in seq_along(percentages$value)){
    if(percentages$value[i]==0){next}
    annot_df = data.frame(
      x = ifelse( percentages[i, g1_score_name], cutoffs[1]*2, cutoffs[1]*0.35) ,
      y = ifelse( percentages[i, g2_score_name], cutoffs[2]*2, cutoffs[2]*0.25) ,
      label = paste0( round(percentages$value[i], 1), "%") )
    if(!is.null(facet_by)) {
      annot_df[[facet_by]] = percentages[i, facet_by]
    }
    p = p + geom_text( data = annot_df, aes(x=x,y=y,label=label) )                
  }
  
  return(p)
}

