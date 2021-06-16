#' Function to read and process taxon metadata file
#'
#' This function reads and process a user specified metadata file or data frame.
#' It assumes that the file is in text format and that the columns are tab-delimited.
#' Metadata provided in any other format other than as a 'character' class file name or
#' a data frame with class "tbl_df", "tbl","grouped_df","data.frame" or "rowwise_df" will not be
#' accepted. When not specified, it will assume that the first column represents
#' the taxon names. The taxon names should match those included in the
#' phylogenetic tree file or object
#'
#' @param taxon.metadata.file Path to the input metadata file name or data frame.
#' @param taxon.metadata.columns A vector containing name of columns in the matadata file or data frame to view in the phylogenetic tree.
#' @param taxon.names A vector containing taxon names to select from the metadata file or data frames. These names must match the taxon names in the phylogenetic tree.
#' @param taxon.id.column Column name in the matadata file or data frame containing taxon names.
#' @param include.first.col A Boolean value specifying whether to use the first column as the taxon names.
#' @param taxon.metadata.delimeter A delimiter separating metadata columns.
#'
#' @return A data frame containing selected metadata columns and taxon names in the phylogenetic tree.
#'
#' @examples
#' \dontrun{
#' Read a tab-delimited file containing metadata
#'
#' metadata.file<-system.file("extdata", "ST320.tsv", package = "RCandy",mustWork = TRUE)
#'
#' metadata.df<-load.taxon.metadata(metadata.file)
#' }
#' @export
#'
#' @import magrittr
#' @import dplyr
#' @importFrom stats setNames
#'
#' @author Chrispin Chaguza, \email{Chrispin.Chaguza@@gmail.com}
#' @references \url{https://github.com/ChrispinChaguza/RCandy}
#'
## Function to read and process taxon metadata file

load.taxon.metadata<-function(taxon.metadata.file,
                              taxon.metadata.columns=NULL,
                              taxon.names=NULL,
                              taxon.id.column=NULL,
                              include.first.col=FALSE,
                              taxon.metadata.delimeter="\t"){
  # Check if the specified metadata file name is a string/character
  if( "character" %in% class(taxon.metadata.file) ){
      # Check if the metadata file exists in the specified path
    if( file.exists(taxon.metadata.file) ){
      # Check if the user has specified the metadata column containing taxon IDs
      tmp.taxon.data<-as_tibble(read.table(taxon.metadata.file,header=TRUE,sep=taxon.metadata.delimeter,comment.char="?",stringsAsFactors=FALSE))
    }else{
      # Stop execution when no valid metadata file is found
      stop("The following metadata file '",taxon.metadata.file,"' was not found.")
    }
  }else{
    if( length(base::setdiff(class(taxon.metadata.file),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
      # Use the data frame name provided
      tmp.taxon.data<-taxon.metadata.file
    }else{
      # Stop execution when no valid metadata file or data frame name is provided
      stop("Something is wrong with the metadata file or data frame. Check that the metadata file
           name is a character or string.")
    }
  }
  # Check that the metadata column names and/or taxon name are provided
  # Convert the taxon names to a named vector if not provided
  if( is.null(taxon.names) ){
    #warning("taxon names are not provided - first column assumed to contain taxon names in the phylogenetic tree.")
    taxon.names<-stats::setNames(1:length(unique(tmp.taxon.data[,1][[1]])),unique(tmp.taxon.data[,1][[1]]))
  }else{
    if( is.null(names(taxon.names)) ){
      taxon.names<-stats::setNames(1:length(unique(taxon.names)),unique(taxon.names))
    }
  }
  # Use all column names in the metadata file or data frame if specific columns to use are not specified
  if( is.null(taxon.metadata.columns) ){
    taxon.metadata.columns<-colnames(tmp.taxon.data)
  }else{
    if( !is.null(taxon.id.column) ){
      # Check if a valid name for the taxon metadata column is specified
      if( !taxon.id.column %in% colnames(tmp.taxon.data) ){
        stop("Column for taxon names not found in the metadata file or data frame")
      }
    }
  }
  # Check if the specified column names are available in the metadata file
  if( length(base::setdiff(taxon.metadata.columns,colnames(tmp.taxon.data)))==0 ){
    if( isTRUE(include.first.col) ){
      if( !is.null(taxon.metadata.columns) ){
        if( is.null(taxon.id.column) ){
          # warning("taxon name column is not specified - first column used.")
          selected.data.columns<-colnames(tmp.taxon.data)
        }else{
          selected.data.columns<-unique(c(taxon.id.column,taxon.metadata.columns))
        }
      }else{
        if( is.null(taxon.id.column) ){
          # warning("No metadata columns specified - all columns used (first column assumed to contain taxon sames)")
          selected.data.columns<-unique(c(colnames(tmp.taxon.data)[1],taxon.metadata.columns))
        }else{
          selected.data.columns<-unique(c(taxon.id.column,taxon.metadata.columns))
        }
      }
      # Select the metadata, assumming the first colum contains taxon IDs
      tmp.data.val<-tmp.taxon.data[tmp.taxon.data[,1][[1]] %in% names(taxon.names),selected.data.columns]
      # Include a column specify the position of each taxon\taxon in the phylogenetic tree
      tmp.data.val$pos<-taxon.names[tmp.data.val[,1][[1]]]
      tmp.data.val<-tmp.data.val[,c("pos",selected.data.columns)]
      return(tmp.data.val)
    }else{
      if( !is.null(taxon.metadata.columns) ){
        if( is.null(taxon.id.column) ){
          # warning("No metadata columns specified - all columns used (first column assumed to contain taxon sames)")
          selected.data.columns<-unique(c(colnames(tmp.taxon.data)[1],taxon.metadata.columns))
        }else{
          selected.data.columns<-unique(c(taxon.id.column,taxon.metadata.columns))
        }
      }else{
        if( is.null(taxon.id.column) ){
          # warning("No metadata columns specified - all columns used (first column assumed to contain taxon sames)")
          selected.data.columns<-colnames(tmp.taxon.data)
        }else{
          selected.data.columns<-unique(c(taxon.id.column,colnames(tmp.taxon.data)))
        }
      }
      # Select the metadata, assumming the first colum contains taxon names
      tmp.data.val<-tmp.taxon.data[tmp.taxon.data[,1][[1]] %in% names(taxon.names),selected.data.columns]
      # Include a column specify the position of each taxon\taxon in the phylogenetic tree
      tmp.data.val$pos<-taxon.names[tmp.data.val[,1][[1]]]
      tmp.data.val<-tmp.data.val[,c("pos",selected.data.columns)]
      return(tmp.data.val)
    }
  }else{
    stop("The following specified metadata columns were not found...",
         paste0(base::setdiff(taxon.metadata.columns,colnames(tmp.taxon.data))))
  }
}

