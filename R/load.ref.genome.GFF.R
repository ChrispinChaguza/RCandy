#' This function loads a reference genome file in GFF format.
#'
#' @param reference.genome Input file name or data frame in GFF file format.
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' Read genome in GFF formatted file, generated usign readseq, and plot the genomic features
#'
#' ref.genome.gff<-system.file("extdata", "Hungary19A-6.gff", package = "RCandy",mustWork = TRUE)
#'
#' new.ref.genome<-load.genome.GFF(ref.genome.gff)
#' }
#'
#' @export
#'
#' @import magrittr
#' @import dplyr
#'
#' @author Chrispin Chaguza, \email{Chrispin.Chaguza@@gmail.com}
#' @references \url{https://github.com/ChrispinChaguza/RCandy}
#'
# Check if a valid reference genome name is provided
load.genome.GFF<-function(reference.genome){
  V1<-seqname<-feature<-score<-strand<-REC<-SEQ<-PROG<-TYPE<-START<-END<-XX<-YY<-ZZ<-NULL
  if( !is.null(reference.genome) & is.character(reference.genome) ){
  # Same coordinates for the genome region to show, default whole genome
  # Read the reference genome GFF annotation file
  tmp.ref.df<-dplyr::as_tibble(read.table(reference.genome,comment.char="#",header=FALSE,sep="\t",fill=TRUE,row.names=NULL,stringsAsFactors = FALSE)) %>%
    dplyr::filter((!grepl("#",V1)) | V1!="seqname" )
  colnames(tmp.ref.df)<-c("seqname","source","feature","start","end","score","strand","frame","attributes")
  tmp.ref.df<-tmp.ref.df[!tmp.ref.df$source %in% c("source"),]
  tmp.ref.df[1,3]<-"source"
  reference.genome.obj1<-tmp.ref.df %>%
    dplyr::mutate(seqname=gsub("# ","",.data$seqname)) %>% mutate(seqname=gsub("^#","",.data$seqname)) %>%
    dplyr::filter(!grepl("gff-version",.data$seqname)) %>%
    dplyr::filter(!.data$feature %in% c("ORIGIN","NA","","##")) %>%
    dplyr::mutate(start=as.integer(.data$start),end=as.integer(.data$end)) %>% dplyr::filter(.data$feature %in% c("source","locus_tag","gene","CDS"))

  # Filter out lines not containing information about the genetic features
  reference.genome.obj<-reference.genome.obj1 %>%
    dplyr::group_by(.data$seqname,.data$source,.data$feature,.data$start,.data$end,.data$score,.data$strand,.data$frame) %>%
    dplyr::filter(!.data$feature %in% c("ORIGIN","NA","","##")) %>%
    dplyr::mutate(attributes=gsub("="," ",.data$attributes)) %>%
    dplyr::mutate(attributes=gsub(":"," ",.data$attributes)) %>% #dplyr::filter(feature %in% c("CDS","gene","source")) %>%
    tidyr::nest() %>% dplyr::mutate(gene=stringr::str_split(stringr::str_trim(stringr::str_split(regmatches(data[[1]],regexpr("(gene|locus_tag|Parent|db_xref|mol_type|organism|ID).*",data[[1]])),";")[[1]][1],side="both")," ")[[1]][-1][1] )
  }else{
    if( length(setdiff(class(reference.genome),c("tbl_df","tbl","data.frame","rowwise_df","grouped_df")))==0 ){
      reference.genome.obj<-reference.genome
    }else{
      # Exit the program when valid genome length is found
      stop("Could not find a feature labelled 'source' in the genome annotation file")
    }
  }
  if( !"source" %in% reference.genome.obj$feature ){
    # Exit the program when valid genome length is found
    stop("Could not find a feature labelled 'source' in the genome annotation file")
  }
  return(reference.genome.obj)
}
