#' Function to read phylogenetic tree in Newick file
#'
#' This function reads and process a user specified phylogenetic tree
#' file in Newick format.
#'
#' @param tree.file.name Path to the input phylogenetic tree file in Newick format.
#'
#' @return A data frame containing selected metadata columns and strain names in the phylogenetic tree.
#'
#' @examples
#' \dontrun{
#' Load phylogenetic tree file.
#'
#' tree.file<-system.file("extdata", "ST320.final_tree.tre", package = "RCandy",mustWork = TRUE)
#'
#' read.tree.file(tree.file.name=tree.file)
#' }
#'
#' @export
#'
#' @import ape
#'
#' @author Chrispin Chaguza, \email{Chrispin.Chaguza@@gmail.com}
#' @references \url{https://github.com/ChrispinChaguza/RCandy}
#'
#### Read and process file or object containing the phylogenetic tree  ####

read.tree.file<-function(tree.file.name){

  # Check if the phylogenetic tree is specified as a phylo object
  if( class(tree.file.name)=="phylo" ){
    tree.to.plot<-tree.file.name

    return(tree.to.plot)
  }else{
    # Check if a valid phylogenetic tree file name is specified
    if( class(tree.file.name)=="character" ){
      # Check if the phylogenetic tree file exists

      if( file.exists(tree.file.name) & !dir.exists(tree.file.name) ){
        tree.to.plot<-read.tree(tree.file.name)

        return(tree.to.plot)
      }else{
        # Stop execution when an invalid phylogenetic tree file or object is specified
        stop("Invalid phylogenetic tree file '",tree.file.name,"' provided...")
      }
    }else{
      # Stop execution when an invalid phylogenetic tree file or object is specified
      stop("Invalid phylogenetic tree file '",tree.file.name,"' provided...")
    }
  }
}
