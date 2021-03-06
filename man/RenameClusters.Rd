% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploratory_pipeline.R
\name{RenameClusters}
\alias{RenameClusters}
\title{Rename clusters using an existing metadata field.}
\usage{
RenameClusters(dge, ident.use = "ident", annot.use = "eday",
  new.name = "ident")
}
\arguments{
\item{dge}{: a Seurat object with field `@scale.data` filled in.}

\item{ident.use}{Cluster labels. The output is the same as this clustering, but renamed.}

\item{annot.use}{embryonic day. For example, if you put eday, a cluster might be called "e11_5_c2".}

\item{new.name}{What field to put the new item in. By default, uses the "ident" slot.}
}
\value{
character vector.
}
\description{
Rename clusters using an existing metadata field.
}
