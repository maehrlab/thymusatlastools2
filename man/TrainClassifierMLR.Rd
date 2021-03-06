% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifier.R
\name{TrainClassifierMLR}
\alias{TrainClassifierMLR}
\title{Train and save a penalized logistic regression classifier.}
\usage{
TrainClassifierMLR(training_dge, results_path, ident.use = "ident",
  genes.use, cells.use = 30, do_knn_smooth = T, k = 30,
  dims.use = 30, reduction.use = "PC", do.save = F, ...)
}
\arguments{
\item{training_dge}{Seurat object for training.}

\item{results_path}{Where to save model.}

\item{ident.use}{Field to use for classifier labels.}

\item{genes.use}{Fields to use as features.}

\item{cells.use}{List of cells or integer N where training is done on N randomly selected cells per training class.}

\item{do_knn_smooth}{Smooth test data via knn before predicting?}

\item{k}{Number of cells to average in knn smoothing}

\item{reduction.use}{@param dims.use For computation of KNN's in a low-dimensional space.}

\item{do.save}{Save model to a file?}

\item{...}{Passed to glmnet.cv.

Uses Seurat::FetchData(training_dge, vars.all = ident.use ) as class labels.
Results (`glmnet` object) and training data (Seurat object) get
saved into a subdirectory of `results_path`. By default, expression values
are smoothed via k-nearest neighbors prior to training, and then training uses 30 cells
randomly selected from each cluster. Sampling is done without replacement unless
fewer than 30 cells are present.}
}
\description{
Train and save a penalized logistic regression classifier.
}
