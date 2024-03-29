% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{SeuratPie}
\alias{SeuratPie}
\title{Make small-multiple pie charts.}
\usage{
SeuratPie(dge, ident.use = "cell_type", facet_by = "eday",
  col = NULL, label = F, main = "Sample makeup by day",
  drop_levels = F)
}
\arguments{
\item{dge}{Seurat object}

\item{ident.use}{Becomes the categories in the pie chart}

\item{facet_by}{Each small multiple contains cases at one level of this variable.}

\item{col}{Optional colorscale.}

\item{label}{Logical. If TRUE, percentages are added.}

\item{main}{Plot title.}

\item{drop_levels}{If TRUE, omit facets that would be empty.}
}
\description{
Make small-multiple pie charts.
}
