% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_expression.R
\name{DiffExpAfterClassifier}
\alias{DiffExpAfterClassifier}
\title{Stratify cells using a classifier, then test for effects of a treatment within each cell type.}
\usage{
DiffExpAfterClassifier(dge_test, dge_train, ident.use, treatment_var,
  treatment_1, treatment_2, ...)
}
\arguments{
\item{dge_test}{Test set for classifier. Seurat object.}

\item{dge_train}{Training set for classifier. Seurat object.}

\item{ident.use}{Metadata field to use as classifier output.}

\item{treatment_var}{Metadata field to use as diff. expression labels.}

\item{treatment_1}{@param treatment_2 Levels of \code{ treatment_var } to contrast.}

\item{...}{Extra args passed to DiffExpByType}
}
\description{
Stratify cells using a classifier, then test for effects of a treatment within each cell type.
}
