#' Check if string is a valid colour name
#'
#' This function checks if a string or a vector contains valid colour names
#'
#' @param colour.name Input character or vector containing colour names
#'
#' @return A Boolean showing whether or not a string or vector contains valid colour names
#'
#' @examples
#' \dontrun{
#' Read the colour names as a string or vector of characters
#'
#' col.names<-is.color(c("red","blue"))
#' col.names<-is.color("blue")
#' }
#'
#' @import grDevices
#'

is.color<-function(colour.name) {

  # Check which vector values are valid colour names
  col.names.check<-sapply(colour.name,function(col.name.val) {
    # Corrently handle exception if invalid color name is encountered
    tryCatch(is.matrix(grDevices::col2rgb(col.name.val)),
             error=function(e) FALSE)
  })

  #Return a named vector showing which items are valid colour names
  return(col.names.check)
}
