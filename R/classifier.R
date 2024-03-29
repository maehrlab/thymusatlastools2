## ------------------------------------------------------------------------

#' Classify cells from one Seurat object in terms of another Seurat object's identity field, with a "reject option" for unfamiliar cells. 
#' 
#' @param dge_train Cells to train classifier on. Seurat object.
#' @param dge_test Cells to be classified. Seurat object.
#' @param ident.use Identity variable to use for training labels.
#' @param vars.all List of raw genes/features to use. If possible, will be accessed 
#' through `FetchData`; in this case, should be numeric. For others, zeroes are filled in.
#' If NULL, uses variable genes from both `dge_train` and `dge_test`.
#' @param my_transform NULL, character, or function. 
#' If `is.null(my_transform)` (default), then `my_transform` is the identity. 
#' if `my_transform` has the form "PCA_<integer>", then the `my_transform` is an unscaled <integer>-dimensional
#' PCA based on the training data. This option triggers special behavior for quantifying 
#' classifier badness, because NN will perform badly in a principal subspace.
#' If a function is given, `my_transform` should accept and return matrices where rows are cells.
#' @param badness Either "pc_dist" or "neighbor_dist" or `NULL`. If `NULL`, default depends on 
#' `my_transform`. You can't use "pc_dist" unless `my_transform` has the form "PCA_<integer>".
#' @param k Number of nearest neighbors to use. Default 25. 
#' @param reject_prop Expected rate of false rejections you're willing to tolerate
#' on held-out training instances. Default is 1/100. This is not honest if `my_transform`
#' is chosen using the training data, and it cannot account for batch effects.
#' @return Seurat object identical to `dge_test` but with new/modified fields for 
#' - `classifier_ident` (predicted class) 
#' - `classifier_badness` (lower means higher confidence)
#' - `classifier_probs_<each identity class from trainset>` (predicted class probabilities)
#' @details Using k-nearest neighbors, classify cells from `dge_test` in terms of 
#' the options in `unique(FetchData(dge_train, ident.use))`, plus a reject option.
#' Rejection happens when the badness (usually distance to the nearest neighbors)
#' falls above a threshold (see `reject_prop`). Badness gets adjusted by cluster,
#' because some clusters naturally are less concentrated on the principal subspace
#' or the coordinates of interest.
#' @export
knn_classifier = function( dge_train, dge_test, ident.use = "ident", 
                           vars.all = NULL, my_transform = "PCA_20", badness = NULL,
                           k = 25, reject_prop = 0 ){
  
  if( is.null( vars.all ) ){ 
    vars.all = union( dge_train@var.genes, dge_test@var.genes ) 
  }
  
  # # Retrieve data, padding test data with zeroes as needed
  coords_train_orig = FetchDataZeroPad( dge_train, vars.all, warn = F )
  coords_test_orig  = FetchDataZeroPad( dge_test,  vars.all, warn = F )
  
  # # set the transform (and badness) 
  if( is.null( my_transform ) ){ my_transform = function(x) x }
  if( is.character( my_transform ) ) { 
    pca_n = strsplit( my_transform, "_", fixed = T )[[1]] 
    atae( "PCA",                  pca_n[1] )
    atat( !is.na( as.numeric( pca_n[2] ) ) )
    cat("Training transform: unscaled PCA of dimension", pca_n[2], " ...\n")
    w = irlba::prcomp_irlba( x = coords_train_orig, n = as.numeric( pca_n[2] ), 
                             retx = F,
                             center = colMeans( coords_train_orig ), 
                             scale = F )
    my_transform = function(x){
      sweep( x, 
             STATS = w$center, 
             FUN = "-",
             MARGIN = 2 ) %*% w$rotation 
    }
    if( is.null( badness ) ){ badness = "pc_dist" }
  } else {
    if( is.null( badness ) ){ badness = "neighbor_dist" }
  }
  if ( badness == "pc_dist" ) {
    cat(" Using distance to principal subspace as badness. \n")
    get_badness = function( nn_dists, x ){ 
      z = sweep( x, 
                 STATS = w$center, 
                 FUN = "-",
                 MARGIN = 2 ) 
      zproj = z %*% w$rotation %*% t( w$rotation )
      square = function( x ) x^2 
      ( z-zproj ) %>% apply( 1, square ) %>% apply( 2, sum ) %>% sapply( sqrt )
    }
  } else {
    cat(" Using average distance to neighbors as badness. \n")
    get_badness = function( nn_dists, x ) { rowMeans( nn_dists )}
  }
  
  cat("Transforming data...\n")
  coords_train_trans = my_transform( as.matrix( coords_train_orig ) ) %>% as.data.frame
  coords_test_trans  = my_transform( as.matrix( coords_test_orig ) ) %>% as.data.frame
  
  # # Find nearest neighbors & classify
  cat("Finding nearest neighbors...\n")
  nn_out = FNN::get.knnx( data = coords_train_trans, 
                          query = coords_test_trans,
                          k=k, algorithm=c( "cover_tree" ) )
  train_labels = FetchData( dge_train, ident.use )[[1]]
  get_label = function(idx) train_labels[idx]
  empir_prob = function(x) factor( x, levels = unique( train_labels ) ) %>% table %>% div_by_sum
  classifier_probs = apply( apply( nn_out$nn.index, 2, get_label ), 1, FUN = empir_prob ) %>% t %>% as.data.frame
  classifier_ident = apply( classifier_probs, 1, function(x) names( which.max( x ) ) )
    
  # # Get badness 
  cat("Calculating badness for each point...\n")
  nn_out_self = FNN::get.knn( data = coords_train_trans,
                              k=k, algorithm=c( "cover_tree" ) )
  classifier_prob_self = apply( apply( nn_out_self$nn.index, 2, get_label ), 
                                1, FUN = empir_prob ) %>% t %>% as.data.frame

  classifier_badness      = get_badness( nn_dists = nn_out$nn.dist,      x = as.matrix( coords_test_orig  ) )
  classifier_badness_self = get_badness( nn_dists = nn_out_self$nn.dist, x = as.matrix( coords_train_orig ) )
  
  # # Adjust badness via regression: some clusters are naturally less dense or farther from PC axes
  model_badness_by_cluster = lm( classifier_badness_self ~ . + 0, 
                                 data = cbind( classifier_badness_self, classifier_prob_self ) )  
  classifier_badness_self = classifier_badness_self - predict( object = model_badness_by_cluster ) 
  classifier_badness      = classifier_badness      - predict( object = model_badness_by_cluster, 
                                                               newdata = classifier_probs ) 
  classifier_badness      = classifier_badness      / sd(classifier_badness_self)
  classifier_badness_self = classifier_badness_self / sd(classifier_badness_self)
  
  # # Set threshold and label rejects
  if( reject_prop > 0 ){ 
    cat("Labeling rejects with attempted controls on false rejection rate ... \n")
    threshold = quantile( classifier_badness_self, 1 - reject_prop )
    hist( classifier_badness_self,  breaks = 80, col = scales::alpha("blue", 0.5))
    hist( classifier_badness,       breaks = 80, col = scales::alpha("red", 0.5),
          add = T, 
          main = "Held-out set badness and threshold",
          xlab = "Average distance to neighbors (train = blue, test = red)" )
    abline( v = threshold )
    text( x = threshold*1.05, y = 100, labels = "Reject", srt = 90 )
    text( x = threshold*0.95, y = 100, labels = "Classify", srt = 90 )
    classifier_ident[classifier_badness >= threshold] = "reject"
  }

  # # Save data and return
  to_add = cbind(          classifier_ident,   
                           classifier_badness,            
                           classifier_probs )
  colnames( to_add ) = c( "classifier_ident", 
                          "classifier_badness" , 
                          "classifier_probs_" %>% paste0( colnames( classifier_probs ) ) )
  rownames( to_add ) = rownames( coords_test_orig )
  dge_test %<>% AddMetaData( to_add )
  
  cat("Done.\n")
  return( dge_test )
}

#' Train and save a penalized logistic regression classifier.
#'
#'
#' @param training_dge Seurat object for training. 
#' @param results_path Where to save model.
#' @param ident.use Field to use for classifier labels. 
#' @param genes.use Fields to use as features.
#' @param cells.use List of cells or integer N where training is done on N randomly selected cells per training class. 
#' @param do_knn_smooth Smooth test data via knn before predicting?
#' @param k Number of cells to average in knn smoothing
#' @param reduction.use @param dims.use For computation of KNN's in a low-dimensional space. 
#' @param do.save Save model to a file?
#' @param ... Passed to glmnet.cv.
#'
#' Uses Seurat::FetchData(training_dge, vars.all = ident.use ) as class labels.
#' Results (`glmnet` object) and training data (Seurat object) get
#' saved into a subdirectory of `results_path`. By default, expression values
#' are smoothed via k-nearest neighbors prior to training, and then training uses 30 cells
#' randomly selected from each cluster. Sampling is done without replacement unless
#' fewer than 30 cells are present.
#' 
#' @export
#' 
TrainClassifierMLR = function( training_dge,
                               results_path, 
                               ident.use = "ident", 
                               genes.use, 
                               cells.use = 30,
                               do_knn_smooth = T, 
                               k = 30,
                               dims.use = 30, 
                               reduction.use = "PC",
                               do.save = F, ... ){

  
  # # Get labels
  training_labels = factor( vectorize_preserving_rownames ( Seurat::FetchData(training_dge, vars.all = ident.use ) ) )

  # # Select cells
  safesample = function(x, n) {
    enough = n <= length(x)
    if( !enough ){
      warning("\nNot enough cells in safesample. Sampling with replacement.\n")
    }
    sample( x, size = n, replace = !enough )
  }
  assertthat::assert_that(!is.factor(cells.use))
  if( is.numeric( cells.use ) ){
    cells.use = aggregate_nice( training_dge@cell.names, training_labels, FUN = safesample, n = cells.use ) %>% c
  }
  
  # # Get features
  genes.use = intersect( genes.use, rownames( training_dge@data ) )
  features_tf = training_dge@data[ genes.use, ] %>% (Matrix::t)

  # # Average the expression values over k nearest neighbors before predicting
  if( do_knn_smooth ){
    cat("\nSmoothing...\n")
    features_smooth = features_tf[cells.use, ]
    embedding = FetchData(training_dge, paste0(reduction.use, 1:dims.use))
    neighbors = FNN::get.knnx( data = embedding, 
                               query = embedding[cells.use, ],
                               k = k, 
                               algorithm = "cover_tree")$nn.index
    for(ii in seq_along( cells.use )){
      features_smooth[ii, ] = Matrix::colMeans(features_tf[neighbors[ii, ], ])
    }
    features_tf = features_smooth
  }
  
  # # Subset labels and imputed expression values
  assertthat::assert_that( all(cells.use %in% names(training_labels)) )
  assertthat::assert_that( all(cells.use %in% rownames(features_tf)) )
  training_labels %<>% extract(cells.use)
  features_tf  %<>% extract(cells.use, )
  
  # # Build classifier (penalized multinomial logistic regression)
  cat("Training classifier...\n")
  # # alpha = 0 does L2 regression. alpha = 1 does LASSO. In between is elastic net.
  mlr_mod = glmnet::cv.glmnet(x = features_tf, y = training_labels, family = "multinomial", ... )
  cat("... classifier trained.\n")
  
  # # Save and return
  if(do.save){
    dir.create.nice( file.path( results_path, "classifier" ) )
    cat("Saving classifier...\n")
    saveRDS(mlr_mod,      file.path( results_path, "classifier", "glmnet_object.Robj" ) )
  }
  cat("\nDone.\n")
  return( mlr_mod )
}

#'
#' @export
#'
train_save_classifier = TrainClassifierMLR

#' Apply a machine learning classifier (from `train_save_classifier`) to new data
#' 
#' @param dge a Seurat object 
#' @param return_type If "seurat" (default), returns a Seurat object with results in the metadata.
#' @param s Controls the regularization used for predictions.
#' @param mlr_mod a glmnet multiclass logistic regression model.
#' You can feed the output of `train_save_classifier` into the `mlr_mod` argument.
#' @param do_knn_smooth Smooth test data via knn before predicting?
#' @param k Number of cells to average in knn smoothing
#' @param reduction.use @param dims.use For computation of KNN's in a low-dimensional space. 
#'
#' @export
#'
ClassifyMLR = function( dge, mlr_mod, return_type = "seurat", s = "lambda.min", 
                        do_knn_smooth = T, k = 30, dims.use = 30, reduction.use = "PC" ){
  # Extract features
  genes.use = ( mlr_mod %>% coef %>% down_idx %>% attributes )$Dimnames %>% down_idx
  features = FetchDataZeroPad( dge, vars.all = setdiff(genes.use, "(Intercept)") ) %>% as.matrix
  
  # Average the expression values over the cell and k-1 nearest neighbors before predicting
  if( do_knn_smooth ){
    features_smooth = features
    neighbors = FNN::get.knn(FetchData(dge, paste0(reduction.use, 1:dims.use)), 
                             k = k - 1, 
                             algorithm = "cover_tree")$nn.index
    for(i in 1:nrow(features) ){
      features_smooth[i, ] = colMeans(features[neighbors[i, ] %>% c(i), ])
    }
    features = features_smooth
  }
  
  # Predict (probabilities)
  predictions = predict( mlr_mod, newx = features, s = s, type = "response" )[, , 1]
  colnames(predictions) %<>% paste0( "classifier_probs_", . )

  # Make hard assignments
  nwm = function(x) names(which.max(x))
  assignments = apply( predictions, 1, nwm )
  assignments %<>% gsub( "classifier_probs_", "", . )
  
  # Assemble data for return
  if(return_type=="seurat"){
    dge %<>% AddMetaData( predictions %>% as.data.frame )
    dge@meta.data$classifier_ident = assignments
    dge %<>% AddClusterIndicators(ident.use = "classifier_ident")
    return( dge )
  } else { 
    return( predictions )
  }
}


#'
#' @export
#'
classify_mlr = ClassifyMLR

#' Summarize results of training a multinomial logistic regression model to classify scRNA data.
#'
#' @param mlr_mod A multinomial logistic regression model from the glmnet package
#' @param dge_train The Seurat object used to train the MLR model. 
#' @param ident.use The field used as the training set labels. 
#' @param results_path Where to put plots.
#' @param plotting_args Plotting params passed to save_feature_plots
#' @param ... Extra params passed to ClassifyMLR
#'
#' @export
#'
DisplayTrainingResults = function( dge_train, results_path, ident.use, mlr_mod, ..., plotting_args = NULL ){
  # Plot cross-validation curve if available
  {
    pdf(file.path(results_path, "glmnet_cv.pdf"), width = 10, height = 7)
    plot(mlr_mod)
    dev.off()
  }
  # Plot all nonzero coefficients by cluster
  ggsave(file.path(results_path, "glmnet_coeffs.pdf"), 
         plot_classifier_coefs( mlr_mod ) + facet_wrap(~variable, ncol = 8, scales = "free"), 
         width = 25, height = 20)

  # ident_levels = levels( dge_train@ident )
  # ident_plus_gene_levels = ident_levels; names(ident_plus_gene_levels) = ident_levels
  # X = get_classifier_coefs( mlr_mod, ... )
  # X$gene = rownames(X)
  # to_plot = reshape2::melt(X, id.vars = "gene") %>% subset(abs(value)>0)
  # for( i in seq_along( ident_levels ) ){
  #   indicator_name = paste0( ident.use, ident_levels[[i]])
  #   label_genes = to_plot %>% subset(select = "gene")
  #   ident_plus_gene_levels[[i]] %<>% paste0( cluster_sizes[[i]], "\n", paste0( label_genes, collapse = "\n") ) 
  # }
  
  # Add annotated levels to Seurat object
  # translate_levels = function(x) ident_plus_gene_levels[as.character( x )]
  # dge@meta.data$ident_plus_markers = translate_levels( ident ) %>% 
  #   factor( ordered = T, levels = translate_levels( ident_levels ) )
    
  # Save coeffs
  write.table(get_classifier_coefs( mlr_mod ), 
              file.path( results_path, "classifier_coeffs.tsv" ),
              sep = "\t", row.names = T, col.names = T, quote =F)
  
  # For measuring classifier performance on training data 
  # df is expected to have original labels and classifier predictions.
  report_concordance = function(df, filename, results_path){
    da = dplyr::arrange
    X = df %>%
      table %>% 
      as.data.frame %>% 
      da(-Freq)
    X$match = as.character(X[[1]])==as.character(X[[2]])
    dir.create.nice(results_path)
    X %>%  write.table(file.path(results_path, filename),
                       quote = F, row.names = F, sep = "\t")
    return(X)
  }
  
  # Classify training set and make sure we are not super underfitted.
  dge_train %<>% ClassifyMLR( mlr_mod, ... )
  dge_train %>% FetchData(c("classifier_ident", ident.use       )) %>% 
    report_concordance( "concordance.txt", results_path = results_path )
  
  # Plot clusters used for training along with genes selected to identify them.
  dge_train %<>% AddClusterIndicators(ident.use = ident.use)
  cc = get_classifier_coefs(mlr_mod)
  for( cluster in colnames(cc) ){
    do.call( what = save_feature_plots, 
             args = list( dge =  dge_train,
                   results_path = results_path, 
                   gene_list =
                     rownames(cc)[which(cc[, cluster]>0)] %>% 
                     c(paste0(ident.use, cluster)) %>%
                     c(paste0("classifier_probs_", cluster)), 
                   gene_list_name = cluster, 
                   types = "PDF" ) %>% c(plotting_args)
    )
  }
}


## ------------------------------------------------------------------------


#' Extract active features for a given class from a glmnet multiclass logistic regression model.
#' 
#' @export
#' 
get_genes_by_cluster = function(mlr_mod, cluster_name, ...){
  X = get_classifier_coefs( mlr_mod, ... )
  assertthat::assert_that(cluster_name %in% colnames( X ) )
  X[cluster_name] %>% vectorize_preserving_rownames %>% is_greater_than(0) %>% which %>% names
} 

#' Extract coefficients from a glmnet multiclass logistic regression model.
#' 
#' @export
#' 
get_classifier_coefs = function( mlr_mod, s = "lambda.min", ... ){
  get_coef = glmnet::coef.cv.glmnet
  cell_types = names( get_coef( mlr_mod, s = s, ... ) )
  genes = rownames( get_coef( mlr_mod, s = s, ... )[[1]] )
  coeffs = data.frame( matrix( NA, nrow = length(genes), ncol = length( cell_types ) ) )
  colnames( coeffs ) = cell_types
  rownames( coeffs ) = make.unique( genes )
  for( cell_type in cell_types ){
    coeffs[[cell_type]] = as.vector( get_coef( mlr_mod, s = s )[[cell_type]] )
  }
  return( coeffs )
}


#' Display coefficients from a glmnet multiclass logistic regression model.
#' 
#' @export
#' 
plot_classifier_coefs = function( mlr_mod, ... ){
  X = get_classifier_coefs( mlr_mod, ... )
  X$gene = rownames(X)
  to_plot = reshape2::melt(X, id.vars = "gene") %>% subset(abs(value)>0)
  to_plot = to_plot[order(to_plot$value), ]
  to_plot$gene_plus_type = interaction(to_plot$gene, to_plot$variable)
  to_plot$gene_plus_type %<>% factor(levels = to_plot$gene_plus_type, ordered = T)
  to_plot$sign = sign(to_plot$value) %>% as.character
  p = ggplot(to_plot) + 
    geom_point( aes( x = value, y = gene_plus_type, colour = sign ) ) + 
    scale_colour_manual(values = c("red", "blue")) +
    facet_wrap(~variable, scales = "free_y") + 
    scale_y_discrete(labels=setNames(to_plot$gene, to_plot$gene_plus_type)) + 
    xlab(label = "weight")
  return(p)
}

#' Given a test set, training set, and classifier, find genes that behave the same in test and train sets.
#'
#' @param dge_test Seurat object
#' @param dge_train Seurat object
#' @param train_ident Metadatum from dge_train used as training labels.
#' @param mlr_mod multiclass logistic regression output from glmnet::cv.glmnet
#' @param ... Additional args passed to get_genes_by_cluster
#'
#' @export
#'
get_active_genes = function( dge_test, dge_train, train_ident, mlr_mod, ... ){
  cell_classes_potential = mlr_mod  %>% get_classifier_coefs(...) %>% colnames
  cell_classes_active    = dge_test %>% FetchData( "classifier_ident" ) %>% extract2(1) %>% unique %>% as.character
  valid_classes = FetchData(dge_train, train_ident ) %>% extract2(1) %>% unique %>% as.character
  if( !all( cell_classes_potential %in% valid_classes ) ){
    stop("Classifier has labels not present in training data. Wrong input for ident.use or dge_train?\n")
  }
  if( !all( cell_classes_active %in% valid_classes ) ){
    stop("Test set has predicted labels not present in training data. Wrong input for dge_test or dge_train?\n")
  }
  positive_genes = lapply( FUN = get_genes_by_cluster, 
                           X = cell_classes_active,
                           mlr_mod = mlr_mod, 
                           ... )
  names(positive_genes) = cell_classes_active
  
  prop_nz_by_cluster = function( dge, genes, ident.use, which.cluster ){
    X = make_heatmap_for_table(dge,
                               desired_cluster_order = FetchData(dge, ident.use )[[1]] %>% unique,
                               genes_in_order = genes, 
                               ident.use = ident.use,
                               return_type = "table", 
                               aggregator = prop_nz,
                               normalize = "none", 
                               norm_fun = function(x) x )
    X = data.frame( gene = rownames(X), 
                    cluster = which.cluster,
                    prop_nz = X[, which.cluster], 
                    stringsAsFactors = F )
    return( X )
  }
  
  prop_nz_by_cluster_vec = function( dge, ident.use ){
    mapply( FUN = prop_nz_by_cluster,
            genes = positive_genes, 
            which.cluster = names( positive_genes ),
            MoreArgs =  list( dge = dge,
                              ident.use = ident.use ), 
            SIMPLIFY = F ) 
  }
  
  positive_genes_expr_test  = prop_nz_by_cluster_vec( dge = dge_test, ident.use = "classifier_ident")
  positive_genes_expr_train = prop_nz_by_cluster_vec( dge = dge_train, ident.use = train_ident )
  merge_one = function( test, train ) { 
    assertthat::are_equal( test[["cluster"]], train[["cluster"]] )
    assertthat::are_equal( test[["gene"]],    train[["gene"]] )
    X = data.frame( gene = test[["gene"]],
                    cluster = test[["cluster"]], 
                    prop_nz_test = test[["prop_nz"]], 
                    prop_nz_train = train[["prop_nz"]], stringsAsFactors = F )
    return(X)
  }
  zipped = mapply( merge_one,
                   test = positive_genes_expr_test, 
                   train = positive_genes_expr_train, 
                   SIMPLIFY = F )
  
  return( zipped )
}



#' Visualize probabilistic classification results
#'
#' @param dge a Seurat object
#' @param results_path folder to place output in
#' @param mlr_mod a glmnet multiclass logistic regression model. Give this or `class_labels`, not both.
#' You can feed the output of `train_save_classifier` into the `mlr_mod` argument.
#' @param class_labels atomic character vector naming columns in the Seurat object that contain
#'  class probabilities. Give this or `mlr_mod`, not both. If names(class_labels) is filled in,
#'  then that's how the spokes get labeled.
#' @param fig_name filename for output.
#' @param facet_by Variable in Seurat object to facet resulting plots by; default is none
#' @param colour_by Variable in Seurat object to map color to; default is none
#' @param fig_name filename for output.
#' @param style If "points", plot each cell as a dot. 
#' If "density", then instead of plotting points, plot 2d density contours.
#' If "hexagons", do AWESOME HEXAGON BINNING YEAHHHHHHH HEXAGONS.
#' @param wheel_order Deprecated.
#' @param do.density Deprecated.
#' @export
wheel_plot = function( dge, results_path, mlr_mod = NULL, class_labels = NULL, fig_name, 
                       colour_by = NULL, facet_by = NULL, style = "points",
                       wheel_order = NULL, do.density = NULL, ... ){
  # # Handle stupid old input
  if(!is.null(wheel_order)){
    warning( "`wheel_order` arg is deprecated. Use `class_labels` instead." )
    if( is.null( class_labels ) ){ class_labels = wheel_order} 
  }
  
  if( !is.null(do.density) && do.density){ 
    warning("do.density is deprecated. Use ` style = \"density\" ` instead.")
    style = "density"
  }
  
  # # Handle regular, non-stupid input
  mlr_mod_given = !is.null( mlr_mod )
  labels_given = !is.null( class_labels )
  if( !mlr_mod_given & !labels_given ){
    stop("Please specify either `mlr_mod` or `class_labels`")
  } else if( !mlr_mod_given & labels_given ) {
    cat( "Fetching stored predictions from Seurat object.\n" )
    predictions = FetchData( dge, vars.all = class_labels ) %>% as.matrix
  } else if( mlr_mod_given & !labels_given ) {
    cat( "Making predictions from glmnet object.\n" )
    predictions = classify_mlr( dge, mlr_mod, return_type = "predictions_only", ... )
    class_labels = colnames( predictions )
  } else if( mlr_mod_given & labels_given ) {
    cat( "Making predictions from glmnet object.\n" )
    warning( "Since `mlr_mod` was given, using `class_labels` 
             only to order wheel spokes, not to fetch stored predictions." )
    predictions = classify_mlr( dge, mlr_mod, return_type = "predictions_only", ... )
    class_labels = colnames( predictions )
  }
  
  # # Make wheel
  lp1 = length( class_labels ) + 1 
  if( is.null( names( class_labels ) ) ){ names( class_labels ) = class_labels }
  unwrapped_circle = seq( 0, 2*pi, length.out = lp1 )
  wheel_spokes = data.frame( z = names( class_labels ), 
                             x = cos( unwrapped_circle[-lp1] ), 
                             y = sin( unwrapped_circle[-lp1] ) )
  # # Process data
  cat("Processing data.\n")
  cell_positions = predictions %*% as.matrix( wheel_spokes[, c("x", "y")] )
  cell_positions = as.data.frame( cell_positions );   colnames( cell_positions ) = c( "x", "y" )
  if( !is.null( colour_by ) ){
      cell_positions[[colour_by]] = Seurat::FetchData(dge, vars.all = colour_by)[ rownames( predictions ),  ] 
  }
  if( !is.null( facet_by ) ){
    cell_positions[[facet_by]] = Seurat::FetchData(dge, vars.all = facet_by )[ rownames( predictions ),  ] 
  }
  
  # # Avoid possible namespace conflict between class_labels in model parameters and class_labels in test data
  if( "class_labels" %in% c( facet_by, colour_by ) ){
    wheel_spokes$class_labels_in_model = wheel_spokes$class_labels
    wheel_spokes$class_labels = NULL
  }
  
  # # Add wheel & label spokes
  wheel_plot = ggplot()   + geom_path ( data = wheel_spokes,   mapping = aes(     x,     y ) )
  wheel_plot = wheel_plot + geom_text ( data = wheel_spokes,   mapping = aes( 1.2*x, 1.2*y, label = z ) )
  # # Add data & facet
  if( style == "density"){
    cat( "Estimating density.\n" )
    wheel_plot = wheel_plot + geom_density_2d( data = cell_positions, 
                                               mapping = aes_string( "x", "y", colour = colour_by ) )
  } else if( style == "hexagons"){
    if( !is.null( colour_by ) ) {
      warning( "Won't use color with style==\"hegaxons\" because fill is mapped to bin count. " )
    }
    wheel_plot = wheel_plot + geom_hex( data = cell_positions, mapping = aes_string( "x", "y" ) ) 
    wheel_plot = wheel_plot +  scale_fill_gradient(trans = "log")
  } else {
    wheel_plot = wheel_plot + geom_point( data = cell_positions, 
                                          mapping = aes_string( "x", "y", colour = colour_by ), alpha = 0.5 ) 
  }
  wheel_plot = wheel_plot + ggtitle( fig_name )
  if( !is.null( colour_by )  && is.numeric( cell_positions[[colour_by]] ) ){
    wheel_plot = wheel_plot + scale_color_gradientn( colours = blue_gray_red )
  }
  if ( !is.null( facet_by ) ){
      wheel_plot = wheel_plot + facet_wrap( as.formula( paste( "~", facet_by ) ) )  
  }
  
  # # Save & return
  cat( "Done. Saving and returning plot.\n" )
  ggsave( filename = file.path( results_path, paste0( fig_name, ".pdf" ) ),
          plot = wheel_plot,
          width = 12, height = 10)
  return( wheel_plot )
}

## ------------------------------------------------------------------------

#' Compute fractional identities for a single cell.
#'
get_frac_ident_pure = function( query_cell, pseudo_bulk ){
  num_pseudo_bulk_datasets = ncol(pseudo_bulk)
  quadprog::solve.QP( Dmat = t(pseudo_bulk) %*% pseudo_bulk,
                      dvec = t(pseudo_bulk) %*% query_cell, 
                      Amat = cbind( rep(1, num_pseudo_bulk_datasets), diag(num_pseudo_bulk_datasets) ),
                      bvec = c(1,                                   rep(0, num_pseudo_bulk_datasets) ), 
                      meq = 1, 
                      factorized = FALSE ) $ solution
}



#' Rank cells by fractional identity and display them with a KDE.
#'
#' @param dge_test @param dge_train Test and training sets (seurat objects).
#' @param main Plot title
#' @param ident.use Field to extract as training labels.
#' @param genes.use Genes to train on.
#' @param test.ident.use @param test.col For plotting the distribution of predictions by group.
#' 
#' @export
#'
display_fractional_identity = function( dge_test, 
                                        dge_train, 
                                        main, 
                                        ident.use, 
                                        genes.use = NULL, 
                                        test.ident.use = NULL, 
                                        test.col = NULL ){
  if( is.null( genes.use ) ){
    dge_train %<>% MeanVarPlot(x.low.cutoff = 0.1, y.cutoff = 0.5)
    genes.use = dge_train@var.genes
  }
  
  # Set up 
  pseudo_bulk = dge_train %>% 
    FetchDataZeroPad(genes.use) %>% 
    aggregate_nice(by = FetchData(dge_train, ident.use )[[1]], FUN = mean) %>% t
  pseudo_bulk
  test_cells = dge_test %>% FetchDataZeroPad(genes.use)
  
  # Obtain the fractional identities via quadratic programming
  fractional_identities = apply( test_cells, 
                                 MARGIN = 1, 
                                 FUN = get_frac_ident_pure,
                                 pseudo_bulk = pseudo_bulk) %>% t %>% as.data.frame
  fi_cols =  colnames( pseudo_bulk )
  colnames( fractional_identities ) = fi_cols
  return( fractional_identities )
}

# #' Given soft assignments between numeric categories, e.g. days of development, rank and plot cells.
# #'
# plot_fractional_identities = function(){
#   test.ident.use = "ident"
#   test.col = DE_IDENT_COLORS
#   fractional_identities[["ranking_var"]] = 
#     as.matrix( fractional_identities[fi_cols] ) %*% (1:length(fi_cols)) %>% 
#     rank( ties.method = "random" )
#   if( !is.null( test.ident.use ) ){
#     fractional_identities[["test_cell_labels"]] = dge_test %>% FetchDataZeroPad( test.ident.use ) %>% extract2(1)
#   }
#   
#   p = ggplot(  fractional_identities[c("ranking_var", "test_cell_labels")] )  + 
#     ggtitle(main) + 
#     geom_density( aes( x = ranking_var, 
#                        colour = make.names(test_cell_labels), 
#                        y = ..scaled.. ) ) +
#     scale_fill_manual  (values = test.col ) +
#     scale_colour_manual(values = test.col ) +
#     xlab( "Cell rank" ) 
#   p = gg_add_colorbar( plot = p, 
#                        x         = fractional_identities[, "ranking_var",      drop = T],
#                        my_labels = fractional_identities[, "test_cell_labels", drop = T] %>% as.character, 
#                        col = test.col )
#   return( plot = p )
# }

