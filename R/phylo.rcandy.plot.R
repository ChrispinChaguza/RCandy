#' Draw and annotate phylogenetic tree with taxon metadata, and
#' recombination events.
#'
#' This function reads a reference genome in GFF format and then plots the
#' genetic features (coding sequences) on both forward and reverse strands.
#'
#' @param tree.file.name File name or ape phylo object containing the phylogenetic tree in Newick format.
#' @param taxon.metadata.file File name or data frame containing taxon metadata.
#' @param taxon.metadata.columns Vector containing metadata columns to plotted.
#' @param gubbins.gff.file Gubbins output recombination file in GFF format.
#' @param recom.input.type Type of input recombination data, either "Gubbins" GFF or "BRATNextGen" tabular data.
#' @param ref.genome.name Reference genome file name in GFF format.
#' @param metadata.column.label.angle Angle for the metadata column labels.
#' @param show.gene.label A Boolean value indicating whether or not to show gene labels in the reference genome.
#' @param gene.label.angle Angle for the gene labels.
#' @param show.metadata.columns A Boolean indicating whether or not to show metadata in the figure.
#' @param subtree.tips A vector containing a subset of taxons/taxa used to generate a subtree from the main phylogenetic tree.
#' @param color.pallette A vector containing names of the viridis colour palletes for visualisation. Choose from "plasma", "cividis", "viridis", "magma" and "inferno".
#' @param taxon.metadata.columns.colors A vector containing column names with custom colours (should be equal in size and same order as the vector specified for taxon.metadata.columns option)
#' @param taxon.id.column Character or string for the column name containing the strain/taxon name in the metadata file.
#' @param show.genome.ticks A Boolean indicating whether to show the xticks for the recombination events diagram/heatmap.
#' @param show.genome.axis A Boolean indicating whether to show the axis for the recombination events diagram/heatmap.
#' @param rec.heatmap.color A two-value vector containing colour names to use for the recombination event diagram/heatmap.
#' @param tree.scale.length A positive number showing the length of the phylogenetic tree branches.
#' @param show.rec.events A Boolean indicating whether to show the recombination event diagram/diagram.
#' @param show.metadata.label A Boolean indicating whether to show labels for the selected metadata columns.
#' @param taxon.metadata.delimeter A delimeter separating metadata columns.
#' @param taxon.metadata.label.cex A number for the size of the labels for the selected matadata columns
#' @param ref.genome.length An optional reference genome length, otherwise it's read from the reference genome GFF file or data frame.
#' @param show.rec.freq.per.base A Boolean indicating whether to show the frequency of recombination per genomic position/base.
#' @param show.rec.freq.per.genome A Boolean indicating whether to show the frequency of recombination events per genome/taxon.
#' @param show.rec.per.genome.scale A boolean indicating whether to show the barplot scale for the number of recombination events per genome.
#' @param rec.events.per.base.as.heatmap A Boolean indicating whether to show the frequency of recombination events per genome/taxon as a barchart or colour scale (heatmap).
#' @param ladderize.tree.right A Boolean indicating whether to ladderize the phylogenetic tree to the right.
#' @param midpoint.root A Boolean indicating whether to root the phylogenetic tree at midpoint.
#' @param rec.plot.bg.transparency A value between 0 and 1 indicating the transparency level of the background for the recombination events diagram/heatmap.
#' @param show.genome.annot A Boolean indicating whether to show genome annotation above the recombination events diagram/heatmap.
#' @param show.rec.plot.tracks A Boolean indicating whether to plot genome tracks for each taxa.
#' @param show.rec.plot.border A Boolean indicating whether to show the border for the recombination events diagram/heatmap.
#' @param ace.model.name A character or string for the model name used for the discrete ancestral character reconstruction. Choose from "ARD", "ER" and "SYM".
#' @param trait.for.ancestral.reconstr A character or string for the column in the metadata file or data frame used for discrete ancestral character reconstruction.
#' @param save.to.this.file If speficified save the plot to this filename, otherwise show the plot in R.
#' @param plot.width Width of the figure
#' @param plot.height Height of the figure
#' @param show.tip.label A Boolean indicating whether to show the phylogenetic tip labels.
#' @param align.tip.label A Boolean indicating whether to align the phylogenetic tip labels.
#' @param show.fig.legend A Boolean indicating whether to show the legend for the selected metadata columns.
#' @param genome.start A positive number indicating start position in the genome to zoom in.
#' @param genome.end A positive number indicating end position in the genome to zoom in.
#' @param color.tree.tips.by.column Character or string for the column name in the metadata file for colouring the phylogenetic tree tips or terminal nodes.
#' @param tree.tip.node.cex A number for the terminal node or tip size in the phylogenetic tree.
#' @param tree.node.cex A number for the size of the nodes phylogenetic tree.
#' @param tree.tip.label.cex A number for the tip label size in the phylogenetic tree.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' Read phylogenetic tree in Newick format, reference genome in GFF formatted file,
#' generated usign readseq, and Gubbins GFF file to plot the genomic features
#'
#' metadata.file<-system.file("extdata", "ST320.tsv", package = "RCandy",
#' mustWork = TRUE)
#' tree.file<-system.file("extdata", "ST320.final_tree.tre", package = "RCandy",
#' mustWork = TRUE)
#' gubbins.gff<-system.file("extdata", "ST320.recombination_predictions.gff",
#' package = "RCandy",mustWork = TRUE)
#' ref.genome.gff<-system.file("extdata", "Hungary19A-6.gff", package = "RCandy",
#' mustWork = TRUE)
#'
#' RCandyVis(tree.file.name = tree.file, taxon.metadata.file = metadata.file,
#' taxon.metadata.columns = c("Source","Country"),ref.genome.name = ref.genome.gff,
#' gubbins.gff.file = gubbins.gff,color.tree.tips.by.column = "Country",
#' show.rec.freq.per.base = FALSE,show.gene.label = FALSE,ladderize.tree.right = TRUE,
#' midpoint.root = TRUE)
#'
#' RCandyVis(tree.file.name = tree.file, taxon.metadata.file = metadata.file,
#' taxon.metadata.columns = c("Source","Country"),ref.genome.name = ref.genome.gff,
#' gubbins.gff.file = gubbins.gff,color.tree.tips.by.column = "Country",
#' show.rec.freq.per.base = FALSE,show.gene.label = FALSE,ladderize.tree.right = TRUE,
#' midpoint.root = TRUE,genome.start = 30000, genome.end = 60000,show.gene.label=TRUE,
#' save.to.this.file = "RCandy.output.pdf",)
#' }
#'
#' @export
#'
#' @import viridis
#' @import ape
#' @import phytools
#' @import shape
#' @import magrittr
#' @import dplyr
#' @import grDevices
#' @import graphics
#' @importFrom stats setNames
#'
#' @author Chrispin Chaguza, \email{Chrispin.Chaguza@@gmail.com}
#' @references \url{https://github.com/ChrispinChaguza/RCandy}
#'
## Generate phylogenetic tree and/or recombination diagram

RCandyVis <- function(tree.file.name,
                      taxon.metadata.file=NULL,
                      taxon.metadata.columns=NULL,
                      gubbins.gff.file=NULL,
                      recom.input.type="Gubbins",
                      ref.genome.name=NULL,
                      metadata.column.label.angle=90,
                      show.gene.label=FALSE,
                      gene.label.angle=45,
                      show.metadata.columns=TRUE,
                      subtree.tips=NULL,
                      color.pallette="inferno",
                      taxon.metadata.columns.colors=NULL,
                      taxon.id.column=NULL,
                      show.genome.ticks=TRUE,
                      show.genome.axis=TRUE,
                      rec.heatmap.color=c("red","blue"),
                      tree.scale.length=NULL,
                      show.rec.events=TRUE,
                      show.metadata.label=TRUE,
                      taxon.metadata.label.cex=0.95,
                      taxon.metadata.delimeter="\t",
                      ref.genome.length=NULL,
                      show.rec.freq.per.base=FALSE,
                      show.rec.freq.per.genome=TRUE,
                      show.rec.per.genome.scale=FALSE,
                      rec.events.per.base.as.heatmap=TRUE,
                      ladderize.tree.right=NULL,
                      midpoint.root=FALSE,
                      rec.plot.bg.transparency=0.10,
                      show.genome.annot=TRUE,
                      show.rec.plot.tracks=FALSE,
                      show.rec.plot.border=FALSE,
                      ace.model.name="ER",
                      trait.for.ancestral.reconstr=NULL,
                      save.to.this.file=NULL,
                      plot.width=12,
                      plot.height=9.5,
                      show.tip.label=FALSE,
                      align.tip.label=FALSE,
                      show.fig.legend=TRUE,
                      genome.start=NULL,
                      genome.end=NULL,
                      color.tree.tips.by.column=NULL,
                      tree.tip.node.cex=0.35,
                      tree.tip.label.cex=0.35,
                      tree.node.cex=0.60){

  # Check if the Boolean arguments are specified correctly
  if(!is.logical(show.gene.label)) stop("'show.gene.label' must be one of TRUE or FALSE")
  if(!is.logical(show.metadata.columns)) stop("'show.metadata.columns' must be one of TRUE or FALSE")
  if(!is.logical(show.genome.ticks)) stop("'show.genome.ticks' must be one of TRUE or FALSE")
  if(!is.logical(show.genome.axis)) stop("'show.genome.axis' must be one of TRUE or FALSE")
  if(!is.logical(show.rec.events)) stop("'show.rec.events' must be one of TRUE or FALSE")
  if(!is.logical(show.metadata.label)) stop("'show.metadata.label' must be one of TRUE or FALSE")
  if(!is.logical(show.rec.freq.per.base)) stop("'show.rec.freq.per.base' must be one of TRUE or FALSE")
  if(!is.logical(show.rec.freq.per.genome)) stop("'show.rec.freq.per.genome' must be one of TRUE or FALSE")
  if(!is.logical(show.rec.per.genome.scale)) stop("'show.rec.per.genome.scale' must be one of TRUE or FALSE")
  if(!is.logical(rec.events.per.base.as.heatmap)) stop("'rec.events.per.base.as.heatmap' must be one of TRUE or FALSE")
  if(!is.numeric(rec.plot.bg.transparency)) stop("'rec.plot.bg.transparency' must be between 0 and 1")
  if(!is.logical(show.genome.annot)) stop("'show.genome.annot' must be one of TRUE or FALSE")
  if(!is.logical(show.rec.plot.border)) stop("'show.rec.plot.border' must be one of TRUE or FALSE")
  if(!is.logical(show.tip.label)) stop("'show.tip.label' must be one of TRUE or FALSE")
  if(!is.logical(align.tip.label)) stop("'align.tip.label' must be one of TRUE or FALSE")
  if(!is.logical(show.fig.legend)) stop("'show.fig.legend' must be one of TRUE or FALSE")

  taxon.metadata.columns.names<-taxon.metadata.columns

  # Check if a valid reference genome name is provided
  if( is.character(ref.genome.name) ){
    # Read the reference genome file
    gb.file1<-load.genome.GFF(ref.genome.name)

    # Extract the reference genome length
    if( "source" %in% gb.file1$feature ){
      ref.genome.length<-gb.file1[gb.file1$feature=="source",]$end
    }else{
      # Stop execution when an invalid reference genome file name is specified
      stop("Could not find a feature labelled 'source' in the genome annotation file")
    }
  }else{
    if( !is.null(ref.genome.name) ){
      if( length(base::setdiff(class(ref.genome.name),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
        if( "source" %in% ref.genome.name$feature ){
          ref.genome.length<-ref.genome.name[ref.genome.name$feature=="source",]$end
        }else{
          # Stop excecution when a valid genome length is found
          stop("Could not find a feature labelled 'source' in the genome annotation file")
        }
      }else{
        # Stop excecution when invalid reference genome file/data is specified
        stop("Something is wrong with the GFF genome annotation file/data/data")
      }
    }
  }

  # Specify the correct start and end of genomic region
  if( is.null(ref.genome.length) ){
    genome.start<-genome.end<-0
  }else{
    if( is.null(genome.start) ){
      genome.start<-0
    }else{
      genome.start<-genome.start
    }
    if( is.null(genome.end) ){
      genome.end<-ref.genome.length
    }else{
      genome.end<-genome.end
    }
    if( genome.start>genome.end | genome.end>ref.genome.length | genome.start>ref.genome.length ){
      # Stop execution when the start genomic position is greater than the end position
      stop("Start position of the selected region should be less than end position in the genome")
    }
    if( genome.start<0 | genome.end<0){
      # Stop execution when the start genomic position is greater than the end position
      stop("Start position of the selected region should not be less than zero")
    }
  }

  # Do not plot the recombination diagram when reference genome length is not specified
  if( is.null(ref.genome.length) ){
    show.genome.annot<-FALSE
    show.rec.events<-FALSE
    genome.start<-genome.end<-NULL
  }else{
    # Check if the start and end position of the genomic region to show is less than the genome length
    if(genome.start>ref.genome.length | genome.end>ref.genome.length){
      # Stop execution otherwise
      stop("Start and end position of the selected region should be less than length of the genome")
    }
  }

  # Check if the genomic region ticks should be shown in the recombination diagram
  if( isTRUE(show.genome.ticks) ){
    show.genome.ticks<-1.0
  }else{
    show.genome.ticks<-0
  }

  # Check if correct input recombination data type is specified
  if( !recom.input.type %in% c("Gubbins","BRATNextGen") ){
    stop("Invalid recombination data specified. Choose from 'Gubbins' or 'BRATNextGen'")
  }

  # Check if the genomic region axis should be shown in the recombination diagram
  if( isTRUE(show.genome.axis) ){
    show.genome.axis<-1.0
  }else{
    show.genome.axis<-0
  }

  # Check the options for plotting the recombination frequency plot
  ##if(!isTRUE(show.rec.freq.per.base) & isTRUE(rec.events.per.base.as.heatmap)){
  ##  warning("'show.rec.freq.per.base' is FALSE...ignoring 'rec.events.per.base.as.heatmap'")
  ##}

  # Specify the viridis colour pallete to use
  if(color.pallette %in% c("viridis","inferno","magma","cividis","plasma")){
    if(color.pallette=="viridis"){
      color.pallette<-viridis::viridis
    }else{
      if(color.pallette=="inferno"){
        color.pallette<-viridis::inferno
      }else{
        if(color.pallette=="magma"){
          color.pallette<-viridis::magma
        }else{
          if(color.pallette=="cividis"){
            color.pallette<-viridis::cividis
          }else{
            if(color.pallette=="plasma"){
              color.pallette<-viridis::plasma
            }
          }
        }
      }
    }
  }else{
    stop("Specified color pallete ",color.pallette," not found. Choose from \"viridis\",\"inferno\",\"magma\",\"cividis\", and \"plasma\"")
  }

  # Check if a valid tree scale length is specified
  if( !is.null(tree.scale.length) ){
    if( (!is.integer(tree.scale.length)) & tree.scale.length<=0 ){
      stop("Phylogenetic tree scale length should be greater than zero...")
    }
  }

  # Read phylogenetic tree object or filename
  # Check if the phylogenetic tree is specified as a phylo object
  if(class(tree.file.name)=="phylo"){
    # Ladderize the phylogenetic tree to the left
    if( is.null(ladderize.tree.right) ){
      # Do not ladderize the phylogenetic tree, leave it as specified
      tree.to.plot<-tree.file.name
    }else{
      # Ladderize the phylogenetic tree to the right
      if( isTRUE(ladderize.tree.right) ){
        tree.to.plot<-ape::ladderize(tree.file.name,right=TRUE)
      }else{
        tree.to.plot<-ape::ladderize(tree.file.name,right=FALSE)
      }
    }
  }else{
    # Ladderize the phylogenetic tree to the left
    if( is.null(ladderize.tree.right) ){
      # Read the tree and leave as it is
      tree.to.plot<-read.tree.file(tree.file.name=tree.file.name)
    }else{
      # Ladderize the phylogenetic tree
      if( isTRUE(ladderize.tree.right) ){
        tree.to.plot<-read.tree.file(tree.file.name=tree.file.name)
        # Midpoint root the phylogenetic tree
        if( isTRUE(midpoint.root) ){
          tree.to.plot<-ape::as.phylo(phytools::midpoint.root(tree.to.plot))
        }
        tree.to.plot<-ape::ladderize(tree.to.plot,right=TRUE)
      }else{
        tree.to.plot<-read.tree.file(tree.file.name=tree.file.name)
        # Midpoint root the phylogenetic tree
        if( isTRUE(midpoint.root) ){
          tree.to.plot<-ape::as.phylo(phytools::midpoint.root(tree.to.plot))
        }
        tree.to.plot<-ape::ladderize(tree.to.plot,right=FALSE)
      }
    }
  }

  # Get the taxon IDs in the same order as specified in the tree

  check.tree.tip <- tree.to.plot$edge[,2] <= length(tree.to.plot$tip.label)
  ordered.tree.tips <- tree.to.plot$edge[check.tree.tip, 2]
  taxons.ordered<-tree.to.plot$tip.label[ordered.tree.tips]

  # Create a vector containing taxon ID and position in the tree
  taxon.names<-stats::setNames(1:length(taxons.ordered),taxons.ordered )

  # Check that metadata file or data frame is available when metadata column for the taxon ID is specified
  if( is.null(taxon.metadata.file) & !is.null(taxon.id.column) ){
    stop("Column for taxa is specified but the metadata file is not provided...")
  }

  # Check that metadata file or data frame is available when metadata column for ancestral reconstruction is specified
  if( is.null(taxon.metadata.file) & !is.null(trait.for.ancestral.reconstr) ){
    stop("Column for taxa is specified but the metadata file is not provided...")
  }

  # Add the specified column for colour metadata among the specified metadata columns
  if( !is.null(color.tree.tips.by.column) ){
    taxon.metadata.columns<-unique(c(taxon.metadata.columns,color.tree.tips.by.column))
  }

  # Check that the metadata file is provided, and warn the use if it's not specified
  if( is.null(taxon.metadata.file) ){
    ##warning("User has not provided metadata file. Phylogenetic tree won't be annotated")
    tmp.data.val<-NULL
  }else{
    # Check if a metadata column for discrete ancestral state reconstruction is specified
    if( !is.null(trait.for.ancestral.reconstr) ){
      # Check if a metadata column for taxon ID is specified
      if( !is.null(taxon.id.column) ){
        # Include the taxon ID column among metadata columns to extract from the metadata file
        taxon.metadata.columns<-unique(c(trait.for.ancestral.reconstr,base::setdiff(taxon.metadata.columns,c(trait.for.ancestral.reconstr,taxon.id.column))))
        # Check if the metadata is specified as a data frame and not a text file
        if( length(base::setdiff(class(taxon.metadata.file),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
          # Check if all the specified metadata columns are available in the data frame
          if( length(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)))==0 ){
            tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                              taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                              taxon.names=taxon.names,taxon.id.column=taxon.id.column,
                                              taxon.metadata.delimeter=taxon.metadata.delimeter)
          }else{
            # Stop execution if some metadata columns are not available in the data frame
            stop("The following columns '",paste(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)),sep="",collapse=", "),"' were not found in the metadata file")
          }
        }else{
          # Read the matadata file and extract the metadata columns
          tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                            taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                            taxon.names=taxon.names,taxon.id.column=taxon.id.column,
                                            taxon.metadata.delimeter=taxon.metadata.delimeter)
        }
      }else{
        # Read metadata file or data frame assumming the taxon ID is in the first column
        ##warning("Column for taxon IDs not specified, assumming it's the first column in the metadata file")
        # Include the taxon ID column among metadata columns to extract from the metadata file
        taxon.metadata.columns<-unique(c(base::setdiff(taxon.metadata.columns,c(trait.for.ancestral.reconstr)),trait.for.ancestral.reconstr))
        if( length(base::setdiff(class(taxon.metadata.file),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
          # Check if some specified metadata columns are not available in the metadata file
          if( length(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)))==0 ){
            tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                              taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                              taxon.names=taxon.names,include.first.col=TRUE,
                                              taxon.metadata.delimeter=taxon.metadata.delimeter)
          }else{
            # Stop execution when some specified metadata columns are not available in the data frame
            stop("The following columns '",paste(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)),sep="",collapse=", "),"' were not found in the metadata file")
          }
        }else{
          # Include the taxon ID column (first column) among metadata columns to extract from the metadata file
          taxon.metadata.columns<-unique(taxon.metadata.columns)
          # Read the metadata file and extract the specified metadata columns
          tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                            taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                            taxon.names=taxon.names,include.first.col=TRUE,
                                            taxon.metadata.delimeter=taxon.metadata.delimeter)
        }
      }
    }else{
      if( !is.null(taxon.id.column) ){
        # Include the taxon ID column (first column) among metadata columns to extract from the metadata file
        taxon.metadata.columns<-c(base::setdiff(taxon.metadata.columns,c(taxon.id.column)))
        if( length(base::setdiff(class(taxon.metadata.file),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
          if(length(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)))==0){
            tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                              taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                              taxon.names=taxon.names,taxon.id.column=taxon.id.column,
                                              taxon.metadata.delimeter=taxon.metadata.delimeter)
          }else{
            # Stop execution when some specified metadata columns are not available in the data frame
            stop("The following columns '",paste(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)),sep="",collapse=", "),"' were not found in the metadata file")
          }
        }else{
          # Read the metadata file and extract the specified metadata columns
          tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                            taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                            taxon.names=taxon.names,taxon.id.column=taxon.id.column,
                                            taxon.metadata.delimeter=taxon.metadata.delimeter)
        }
      }else{
        ##warning("Column for taxon IDs not specified, assumming it's the first column in the metadata file")
        if( length(base::setdiff(class(taxon.metadata.file), c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
          if(length(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)))==0){
            # Include the taxon ID column (first column) among metadata columns to extract from the metadata file
            taxon.metadata.columns<-unique(c(taxon.metadata.columns))
            tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                              taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                              taxon.names=taxon.names,include.first.col=TRUE,
                                              taxon.metadata.delimeter=taxon.metadata.delimeter)
          }else{
            # Stop execution when some specified metadata columns are not available in the data frame
            stop("The following columns '",paste(base::setdiff(c(taxon.metadata.columns),colnames(taxon.metadata.file)),sep="",collapse=", "),"' were not found in the metadata file")
          }
        }else{
          # Read the metadata file and extract the specified metadata columns
          tmp.data.val<-load.taxon.metadata(taxon.metadata.file=taxon.metadata.file,
                                            taxon.metadata.columns=unique(c(taxon.metadata.columns,taxon.metadata.columns.colors)),
                                            taxon.names=taxon.names,include.first.col=TRUE,
                                            taxon.metadata.delimeter=taxon.metadata.delimeter)
        }
      }
    }
  }

  # Check that the number of taxon/taxa in the metadata data frame is the same as in the phylogenetic tree
  if( !is.null(tmp.data.val) & !is.null(taxon.metadata.file) ){
    if( length(base::setdiff(names(taxon.names),c(unname(unlist(tmp.data.val[,2]))) ))!=0 ){
      stop("Some taxon names were not found in the metadata file or data frame")
    }
  }
  # Hide the colour strips and figure legend if metadata columns are not specified
  if( is.null(taxon.metadata.columns) | is.null(tmp.data.val) | is.null(taxon.metadata.file)  ){
    show.metadata.columns<-FALSE
    show.fig.legend<-FALSE
  }
  # Check if Gubbins recombination file is provided and valid
  if( !is.null(gubbins.gff.file) ){
    if( is.character(gubbins.gff.file) ){
      if( file.exists(gubbins.gff.file) ){
        tree.rec.data<-load.gubbins.GFF(gubbins.gff.file,recom.input.type=recom.input.type)
      }else{
        stop("Cannot find Gubbins file '",gubbins.gff.file,"' containing recombination events")
      }
    }else{
      if( length(base::setdiff(class(gubbins.gff.file),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0 ){
        tree.rec.data<-gubbins.gff.file
      }else{
        stop("Incorrect format for the Gubbins file '",gubbins.gff.file,"' is specified...")
      }
    }
  }else{
    ##warning("Gubbins file containing recombination events not provided...")
    tree.rec.data<-NULL
  }

  # Check if the taxon names for the subtree are specified correctly
  if( !is.null(subtree.tips)){
    if( !is.vector(subtree.tips) ){
      stop("Argument for 'subtree.tips' should be a vector")
    }
  }

  # Check if specified tips for subtree are available in the full tree
  if( !is.null(subtree.tips) ){
    # Check if the correct number of taxons/taxa for the subtree is specified
    if( length(subtree.tips)<2 ){
      stop("Subtree should contain >=2 taxons/taxa")
    }
    # Extract subtree from the specified taxon IDs
    tree.to.plotr<-drop.tip( tree.to.plot,tip=base::setdiff( tree.to.plot$tip.label,c( subtree.tips ))    )
    # Extract recombinatione events involving the taxons/taxa specified in the subtree
    tmp.rec.df<-tree.rec.data[1,][-1,]
    for(i in 1:length(tree.rec.data$SEQ)){
      tmp.gene.val<-c()
      for(j in unlist(tree.rec.data$gene[i]) ){
        if( j %in% tree.to.plotr$tip.label ){
          tmp.gene.val<-c(tmp.gene.val,j)
        }
      }
      tmp.gene.val<-unique(tmp.gene.val)
      if( is.null(tmp.rec.df) ){
        if( length(tmp.gene.val)!=0 ){
          tmp.rf<-tree.rec.data[i,]
          tmp.rf$gene<-list(tmp.gene.val)
          tmp.rec.df<-tmp.rf
        }
      }else{
        if( length(tmp.gene.val)!=0 ){
          tmp.rf<-tree.rec.data[i,]
          tmp.rf$gene<-list(tmp.gene.val)
          tmp.rec.df<-rbind(tmp.rec.df,tmp.rf)
        }
      }
    }
    # Check that the specified taxon IDs are available in the full tree
    if( length(base::setdiff(subtree.tips,tree.to.plotr$tip.label))!=0 ){
      stop("Some selected taxa for subtree not found in the full tree")
    }
    # Update the variables containg the phylogenetic tree and recombination data to plot
    tree.to.plot<-tree.to.plotr
    tree.rec.data<-tmp.rec.df
    # Get ordered list of taxon IDs from the new tree (subtree)
    check.tree.tip <- tree.to.plot$edge[,2] <= length(tree.to.plot$tip.label)
    ordered.tree.tips <- tree.to.plot$edge[check.tree.tip, 2]
    taxons.ordered<-tree.to.plot$tip.label[ordered.tree.tips]
    taxon.names<-stats::setNames(1:length(taxons.ordered),taxons.ordered )
    tmp.data.val<-tmp.data.val[unname(unlist(tmp.data.val[,2])) %in% taxons.ordered,] %>% as.data.frame()
    temp.position<-stats::setNames(1:length(taxons.ordered),taxons.ordered)
    tmp.data.val$pos<-unname(unlist(temp.position[ unname(unlist(tmp.data.val[,2])) ]))
    rownames(tmp.data.val)<-unname(unlist(tmp.data.val[,2]))
    tmp.data.val<-tmp.data.val[taxons.ordered,]
  }

  # Check if the user specified colors are valid
  if(!is.null(taxon.metadata.columns.colors)){
    if(is.null(color.tree.tips.by.column)){
      if(length(taxon.metadata.columns.names)!=length(taxon.metadata.columns.colors)){
        stop("Number of metadata columns and columns containing custom colour definitions are not equal")
      }
    }else{
      taxon.metadata.columns.id<-setdiff(taxon.metadata.columns,color.tree.tips.by.column)
      if((length(taxon.metadata.columns.names)!=length(taxon.metadata.columns.colors))){
        stop("Number of metadata columns and columns containing custom colour definitions are not equal")
      }
    }
  }

  # Check that valid custom colours are provided for the recombination diagram
  if( length(rec.heatmap.color)==2 ){
    if( !is.color(rec.heatmap.color[1]) & is.color(rec.heatmap.color[2]) ){
      stop("Invalid heatmap colour name specified...")
    }
  }else{
    stop("Invalid number of heatmap colours specified...")
  }

  # Set default values for different panels in the layout canvas
  strip.legend.size<-0.30
  metadata.panel.width<-0.15
  annot.label.height<-0.60
  annot.panel.height<-0.20
  rec.events.panel<-1.00
  tree.panel.width<-2.50
  last.panel.width<-0.30
  first.panel.width<-0.15
  tree.panel.height<-2.50
  bottom.panel.width<-0.15

  # Specify the correct dimensions for the genome annotation panels
  if( isTRUE(show.genome.annot) & !is.null(gubbins.gff.file) ){
    # Specify correct dimensions for the genome annotation panels when recombination and genome annotations are available
    if(  isTRUE(show.gene.label) & isTRUE(show.genome.annot) ){
      annot.panel.height<-0.40
      annot.label.height<-0.60
    }else{
      if( isTRUE(show.metadata.label) ){
        annot.label.height<-0.30
      }else{
        annot.label.height<-0.30
      }
      annot.panel.height<-0.30
    }
  }else{
    # Specify correct dimensions for the genome annotation panels when either recombination or genome annotations are not available
    if( isTRUE(show.gene.label) & isTRUE(show.genome.annot) ){
      annot.panel.height<-0.40
      annot.label.height<-0.60
    }else{
      if( isTRUE(show.metadata.label) ){
        annot.label.height<-0.40
      }else{
        annot.label.height<-0.40
      }
      annot.panel.height<-0.40
    }
  }

  # Specify correct dimension for the recombination heatmap
  if( !is.null(ref.genome.name) & !is.null(gubbins.gff.file) & isTRUE(show.rec.events) ){
    rec.events.panel<-3
  }else{
    rec.events.panel<-0.05
  }

  # Check if the plot should be displayed within R or saved to a file
  if( !is.null(save.to.this.file) ){
    pdf(file=save.to.this.file,width=plot.width,height=plot.height)
  }

  taxon.metadata.columns.id<-c(taxon.metadata.columns.names,
                               setdiff(c(color.tree.tips.by.column,trait.for.ancestral.reconstr),taxon.metadata.columns.names))

  # Check if metadata colour strips and figure legend should be plotted and specify the correct dimension for the panels
  if( !is.null(taxon.metadata.file) ){
    if(length(taxon.metadata.columns.names)<=4 ){
      strip.legend.size<-0.50
      metadata.panel.width<-0.50
    }else{
      if( isTRUE(show.metadata.columns.names) ){
        metadata.panel.width<-(length(taxon.metadata.columns.names)/4)*1.05
      }else{
        metadata.panel.width<-0.05
      }
    }
  }else{
    metadata.panel.width<-0.05
    strip.legend.size<-0.05
  }
  if( isTRUE(show.fig.legend) ){
    strip.legend.size<-0.50+((length(taxon.metadata.columns.id)/4)*0.20)
  }else{
    strip.legend.size<-0.20
  }

  # Specify correct panel for the recombination frequency/heatmap per genomic position
  if( isTRUE(show.rec.freq.per.base) ){
    if(!isTRUE(rec.events.per.base.as.heatmap)){
      rec.heatmap.size<-0.25
    }else{
      rec.heatmap.size<-0.125
    }
  }else{
    rec.heatmap.size<-0.25
  }

  # Specify correct panel for the recombination frequency/heatmap per genome
  if( isTRUE(show.rec.freq.per.genome) ){
    rec.freq.per.taxon.panel<-0.50
  }else{
    rec.freq.per.taxon.panel<-0.05
  }

  # Specify the correct dimensions for the panel between metadata columns and recombination diagram
  space.panel.width<-0.05

  # Generate layout of all the plots based on the specified dimensions
  {
    layout(mat=matrix(c(1,2,3,13,10,19,16,8,
                        1,2,3,13,7,19,17,8,
                        1,4,5,13,6,19,15,8,
                        1,9,9,13,11,19,18,8,
                        1,12,12,12,12,12,12,8,
                        14,14,14,14,14,14,14,8),nrow=6,byrow=TRUE),
           heights=c(annot.label.height,annot.panel.height,tree.panel.height,rec.heatmap.size,strip.legend.size,bottom.panel.width),
           widths=c(first.panel.width,tree.panel.width,metadata.panel.width,space.panel.width,rec.events.panel,space.panel.width,rec.freq.per.taxon.panel,last.panel.width))
  }

  ################################  PANEL 1 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 2 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 3 - COLOR COLUMN STRIP NAMES  #####################################
  {
    ## Show metadata column colour strips in the diagram
    # Check if metadata column colour strip labels should be shown and whether the metadata file or data frame is available

    taxon.metadata.columns.id<-taxon.metadata.columns.names

    if( isTRUE(show.metadata.columns) & !is.null(taxon.metadata.file) & isTRUE( show.metadata.label ) ){
      # Add metadata column names after the phylogenetic tree
      par(mai=c(0,0,0,0))
      plot(1:4,1:4,col=rgb(0,0,0,alpha=0),bty="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA,
           xlim=c(0,length(taxon.metadata.columns.id)+0.5),ylim=c(1,length(taxons.ordered)),
           xaxs="i",yaxs="r")
      loop.val<-1
      # Add each metadata column name
      for(count.val in taxon.metadata.columns.id){
        text(loop.val,1,count.val,srt=metadata.column.label.angle,adj=0,cex=taxon.metadata.label.cex)
        loop.val<-loop.val+1
      }
    }else{
      # Otherwise hide metadata column names after the phylogenetic tree
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }

  ################################  PANEL 4 - PHYLOGENETIC TREE  #####################################
  {
    ## Show phylogenetic tree in the diagram
    # Check if metadata column for ancestral reconstruction is specified
    if( !is.null(trait.for.ancestral.reconstr) ){
      # Perform discrete ancestral character reconstruction using the specified metadata column
      if( (ace.model.name %in% c("ARD","ER","SYM")) & (length(ace.model.name)==1) ){
        # Halt execution if the metadata file or data frame is not provided
        if( is.null(taxon.metadata.file) | is.null(tmp.data.val) ){
          stop("Metadata file '",taxon.metadata.file,"' not found...")
        }else{
          # Run ancestral character reconstruction using the default or specified model
          # Transform zero tree branch lengths to very small values to avoid errors
          tree.to.plot$edge.length[tree.to.plot$edge.length==0]<-max(nodeHeights(tree.to.plot))*1e-6
          # Prepare the data for acnestral reconstruction
          tmp.data.val.X<-tmp.data.val[tmp.data.val[,2][[1]] %in% tree.to.plot$tip.label,]
          pheno.to.reconstr<-stats::setNames(unname(unlist(tmp.data.val.X[,trait.for.ancestral.reconstr])),unname(unlist(tmp.data.val.X[,2])) )
          pheno.to.reconstr<-as.factor(pheno.to.reconstr)
          cols<-stats::setNames(color.pallette(length(unique(pheno.to.reconstr))),levels(pheno.to.reconstr))

          # Run ancestral reconstruction
          fitARD.ALL<-NULL
          if( length(unique(pheno.to.reconstr))>1 ){
            ancestral.reconstr.warning <- tryCatch(fitARD.ALL<-ace(x=pheno.to.reconstr,phy=tree.to.plot,model=ace.model.name,type="discrete",CI=TRUE),
                                                   error=function(e) e, warning=function(w) w)
            if(methods::is(ancestral.reconstr.warning,"warning")){
              stop("Too many states may have been specified for ancestral reconstruction but with insufficient tip/taxon data. Try specifying a different model for the 'ace.model.name' parameter")
              }
          }else{
            trait.for.ancestral.reconstr<-NULL
            cat("Skipping ancestral character reconstruction on the phylogeny - at least two trait values are required")
          }
        }
      }else{
        # Halt execution if an invalid model is specified
        stop("Model for discrete ancestral reconstruction not found or correctly specified...")
      }
    }

    # Plot phylogenetic tree
    if( isTRUE(show.metadata.columns) & !is.null(taxon.metadata.file) ){
      par(mai=c(0,0,0,0))
      if( isTRUE(align.tip.label) ){
        # Align tip labels in the tree
        if( is.null(trait.for.ancestral.reconstr) ){
          # Do not include ancestral character reconstruction data
          plot.phylo(tree.to.plot,show.tip.label=TRUE,align.tip.label=align.tip.label,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
        }else{
          plot.phylo(tree.to.plot,show.tip.label=TRUE,align.tip.label=align.tip.label,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
          if( length(unique(pheno.to.reconstr))>1 ){
            # Include ancestral character reconstruction data
            ape::nodelabels(pie=fitARD.ALL$lik.anc,cex=tree.node.cex,piecol = cols)
            ape::tiplabels(pie=to.matrix(pheno.to.reconstr[tree.to.plot$tip.label],levels(pheno.to.reconstr)),piecol=cols,cex=tree.tip.node.cex)
          }
        }
      }else{
        # Do not align tip labels in the tree
        if( is.null(trait.for.ancestral.reconstr) ){
          # Do not include ancestral character reconstruction data
          plot.phylo(tree.to.plot,show.tip.label=show.tip.label,align.tip.label=FALSE,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
        }else{
          plot.phylo(tree.to.plot,show.tip.label=show.tip.label,align.tip.label=FALSE,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
          if( !is.null(taxon.metadata.file) & length(unique(pheno.to.reconstr))>1 ){
            if( length(unique(pheno.to.reconstr))>1 ){
              # Include ancestral character reconstruction data
              ape::nodelabels(pie=fitARD.ALL$lik.anc,cex=tree.node.cex,piecol = cols)
              ape::tiplabels(pie=to.matrix(pheno.to.reconstr[tree.to.plot$tip.label],levels(pheno.to.reconstr)),piecol=cols,cex=tree.tip.node.cex)

            }
          }
        }
      }
    }else{
      par(mai=c(0,0,0,0))
      if( isTRUE(align.tip.label) ){
        # Align tip labels in the tree
        if( is.null(trait.for.ancestral.reconstr) ){
          # Do not include ancestral character reconstruction data
          plot.phylo(tree.to.plot,show.tip.label=TRUE,align.tip.label=align.tip.label,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
        }else{
          plot.phylo(tree.to.plot,show.tip.label=TRUE,align.tip.label=align.tip.label,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
          if( !is.null(taxon.metadata.file)  ){
            # Include ancestral character reconstruction data
            ape::nodelabels(pie=fitARD.ALL$lik.anc,cex=tree.node.cex,piecol = cols)
            ape::tiplabels(pie=to.matrix(pheno.to.reconstr[tree.to.plot$tip.label],levels(pheno.to.reconstr)),piecol=cols,cex=tree.tip.node.cex)
          }
        }
      }else{
        if( is.null(trait.for.ancestral.reconstr) ){
          # Do not include ancestral character reconstruction data
          plot.phylo(tree.to.plot,show.tip.label=show.tip.label,align.tip.label=FALSE,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
        }else{
          plot.phylo(tree.to.plot,show.tip.label=show.tip.label,align.tip.label=FALSE,
                     cex=tree.tip.label.cex,edge.width=1.251,yaxs="r")
          if( !is.null(taxon.metadata.file) ){
            # Include ancestral character reconstruction data
            ape::nodelabels(pie=fitARD.ALL$lik.anc,cex=tree.node.cex,piecol = cols)
            ape::tiplabels(pie=to.matrix(pheno.to.reconstr[tree.to.plot$tip.label],levels(pheno.to.reconstr)),piecol=cols,cex=tree.tip.node.cex)
          }
        }
      }
    }
    if( !is.null(color.tree.tips.by.column) ){
      if(!is.null(color.tree.tips.by.column) & !is.null(trait.for.ancestral.reconstr)){
        cat("Tips will be coloured by trait.for.ancestral.reconstr option and color.tree.tips.by.column will be ignored")
      }else{
        if( color.tree.tips.by.column %in% taxon.metadata.columns ){
          tmp.data.val.X1<-tmp.data.val[unname(unlist(tmp.data.val[,2])) %in% tree.to.plot$tip.label,]
          tmp.data.pos<-unname(unlist(tmp.data.val.X1[,1]))
          tmp.data.trait<-tmp.data.val.X1[,color.tree.tips.by.column]
          tmp.data.id<-unname(unlist(tmp.data.val.X1[,2]))
          tmp.dat.final<-data.frame(rank=tmp.data.pos,id=tmp.data.id,trait=tmp.data.trait)
          rownames(tmp.dat.final)<-tmp.dat.final$id
          tmp.dat.final<-tmp.dat.final[tree.to.plot$tip.label,]
          tmp.dat.final$rank<-1:length(tmp.dat.final$rank)
          tree.tip.vals<-tmp.dat.final[,c(1,2,3)] %>% dplyr::mutate(VAL=1) %>% dplyr::arrange(3) %>%
            tidyr::spread(key=3,value=4,fill=0) %>% dplyr::select(-id) %>% dplyr::arrange(.data$rank) %>%
            tibble::column_to_rownames(var="rank") %>% as.matrix()
          tree.tip.vals<-tree.tip.vals/rowSums(tree.tip.vals)
          ape::tiplabels(tip=1:dim(tree.tip.vals)[1],pie=tree.tip.vals,
                         piecol=color.pallette(dim(tree.tip.vals)[2]),cex=tree.tip.node.cex)
        }else{
          stop("Attribute for colouring tree tips not found in the metadata")
        }
      }
    }
  }

  ################################  PANEL 5 - COLUMN COLOUR STRIPS  #####################################
  {
    ## Show metadata column colour strips
    # Check if the metadata column colour strips should be shown and whether the metadata file or data frame is available
    if( isTRUE(show.metadata.columns) & !is.null(taxon.metadata.file) ){
      # Add metadata column colour strips after the phylogenetic tree

      taxon.metadata.columns.id<-taxon.metadata.columns.names

      par(mai=c(0,0,0,0))
      plot(1,0.5,col=rgb(0,0,0,alpha=0),bty="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA,
           xlim=c(0,length(taxon.metadata.columns.id)+0.5),ylim=c(1,length(taxons.ordered)+0.5),xaxs="i",yaxs="r" )
      loop.val<-1

      for(count.val in taxon.metadata.columns.id){
        if(is.null(taxon.metadata.columns.colors)){
          strips.tmp<-tmp.data.val[order(tmp.data.val$pos),c("pos",count.val)]
          strips.vals<-gsub("^NA$","N/A",sort(unique(unname(unlist(strips.tmp[,count.val])))))
          strips.tmp.cols<-stats::setNames(color.pallette(length(strips.vals)),strips.vals )
          strips.tmp$col<-strips.tmp.cols[ sapply(unname(unlist(strips.tmp[,2])), as.character)  ]
          strips.tmp<-strips.tmp %>% unique()
          colnames(strips.tmp)<-c("pos","trait","col")
        }else{
          strips.tmp<-data.frame(col=tmp.data.val[order(tmp.data.val$pos),taxon.metadata.columns.colors[loop.val]][[1]],
                                 pos=tmp.data.val[order(tmp.data.val$pos),"pos"][[1]]) %>%
            dplyr::rowwise() %>% dplyr::mutate(col=ifelse(isTRUE(unname(is.color(.data$col))),col,NA))
        }

        graphics::rect(loop.val-0.5,strips.tmp$pos-0.5,
                       loop.val+0.5,strips.tmp$pos+0.5,col=strips.tmp$col,lwd=0.001,border=strips.tmp$col)
        loop.val<-loop.val+1
      }
    }else{
      # Hide the metadata column colour strips after the phylogenetic tree
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }

  ################################  PANEL 6 - RECOMBINATION EVENTS PLOT #####################################
  {
    ## Plot location of recombination events in the genome
    # Check if the data frame containing recombination events is available
    if( !is.null(gubbins.gff.file) ){
      if( isTRUE(show.rec.events) ){
        # Show recombination events in the diagram
        par(mai=c(0,0,0,0))
        plot(0,0,las=1,xlim=c(genome.start,genome.end),ylim=c(1,length(taxon.names)+0.5 ),bty="n",xaxt="n",yaxt="n",
             xlab="",ylab="Genome",col=rgb(0,0,0,alpha=0),xaxs="i",yaxs="r")
        # Show x-axis only if recombination frequency plot is not shown

        genome.ticks <- pretty(c(genome.start,genome.end))

        ##genome.ticks <- formatC(seq(genome.start,genome.end,floor((genome.end-genome.start+1)/5)), format="d", digits=0)
        if( !isTRUE(show.rec.freq.per.base) ){
          if( isTRUE(show.rec.freq.per.genome) ){
            axis(1,at=genome.ticks,lwd.ticks=show.genome.ticks,lwd=show.genome.axis,
                 labels=formatC(pretty(c(genome.start, genome.end)),digits = 0,format = "d" ))
          }else{
            axis(1,at=genome.ticks,lwd.ticks=show.genome.ticks,lwd=show.genome.axis,
                 labels=formatC(pretty(c(genome.start, genome.end)),digits = 0,format = "d" ))
          }
        }
        # Show background colour for the recombination matrix diagram
        if( is.numeric(rec.plot.bg.transparency) & !is.null(ref.genome.name) ){
          graphics::rect(genome.start-0.5,0.45,
                         genome.end,length(taxon.names)+0.45,lty=6,col=rgb(0,0,0,alpha=rec.plot.bg.transparency),lwd=0.0,border=rgb(0,0,0,alpha=0.035))
        }
        # Adjust the line colours appropriately depending on whether the diagram is shown in R or saved to file
        if( is.null(save.to.this.file) ){
          for(count.val in 1:length(tree.rec.data$SEQ)){
            graphics::rect(tree.rec.data[count.val,]$START,taxon.names[ tree.rec.data[count.val,]$gene[[1]]]-0.50,
                           tree.rec.data[count.val,]$END,taxon.names[ tree.rec.data[count.val,]$gene[[1]]]+0.50,
                           border=rgb(1,1,1),lwd=0.001,
                           col=ifelse(length(taxon.names[ tree.rec.data[count.val,]$gene[[1]] ])>1,rec.heatmap.color[1],rec.heatmap.color[2]))
          }
        }else{
          for(count.val in 1:length(tree.rec.data$SEQ)){
            graphics::rect(tree.rec.data[count.val,]$START,taxon.names[ tree.rec.data[count.val,]$gene[[1]]]-0.50,
                           tree.rec.data[count.val,]$END,taxon.names[ tree.rec.data[count.val,]$gene[[1]]]+0.50,
                           border=rgb(1,1,1),lwd=0.001,
                           col=ifelse(length(taxon.names[ tree.rec.data[count.val,]$gene[[1]] ])>1,rec.heatmap.color[1],rec.heatmap.color[2]))
          }
        }
        # Check if genome tracks should be shown for each taxon/taxa
        if( isTRUE(show.rec.plot.tracks) ){
          graphics::rect(genome.start-0.5,(0:length(taxon.names))-0.05-0.5,
                         genome.end,(0:length(taxon.names))+0.05-0.5,
                         lwd=1,col=rgb(1,1,1),border=rgb(1,1,1))
        }
        if( isTRUE(show.rec.plot.border) ){
          graphics::rect(genome.start-0.5,0.5,
                         genome.end,length(taxon.names)+0.5,
                         col=rgb(0,0,0,alpha=0.035),lwd=1.0,border=rgb(0,0,0,alpha=1.0))
        }
      }else{
        # Do not show recombination events in the diagram
        # Empty plot as a placeholder
        par(mai=c(0,0,0,0))
        show.blank.plot()
      }
    }else{
      # Do not show recombination events in the diagram
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }

  ################################  PANEL 7 - GENE FEATURES  #####################################
  {
    ## Plot genome annotation features
    # Check if recombination events and genome annotation data is available
    if( !is.null(gubbins.gff.file) & !is.null(ref.genome.name) ){
      # Check if the genome annotations should be plotted
      if( isTRUE(show.rec.events) ){
        # Show genome annotations
        if( isTRUE(show.genome.annot) & !is.null(taxon.metadata.file) ){
          par(mai=c(0,0,0,0))
          if( is.null(save.to.this.file) ){
            show.genome.annotation.plot(genome.name=ref.genome.name,
                                        genome.start=genome.start,genome.end=genome.end,
                                        show.gene.label=FALSE,genome.start.upstream=0,genome.end.downstream=0,gene.feature.width=1.05)
          }else{
            show.genome.annotation.plot(genome.name=ref.genome.name,
                                        genome.start=genome.start,genome.end=genome.end,
                                        show.gene.label=FALSE,genome.start.upstream=0,genome.end.downstream=0,gene.feature.width=0.85)
          }
        }else{
          # Hide genome annotations
          # Empty plot as a placeholder
          par(mai=c(0,0,0,0))
          show.blank.plot()
        }
      }else{
        # Hide genome annotations
        # Empty plot as a placeholder
        par(mai=c(0,0,0,0))
        show.blank.plot()
      }
    }else{
      # Hide genome annotations
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }

  ################################  PANEL 8 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 9 - TREE SCALE BAR #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    plot.phylo(tree.to.plot,show.tip.label=FALSE,align.tip.label=FALSE,
               cex=0.55,edge.width=0.0001,yaxs="r",edge.color = rgb(1,1,1,alpha=0.0))
    if( !is.null(tree.scale.length) ){
      ape::add.scale.bar(0,40,cex=1,col="black",lwd=1.5,length=tree.scale.length)
    }else{
      ape::add.scale.bar(0,40,cex=1,col="black",lwd=1.5)
    }
  }

  ################################  PANEL 10 - GENE FEATURES LABELS  #####################################
  {
    ## Show genome annotation feature labels
    # Check if recombination events are shown in the diagram
    if( isTRUE(show.rec.events) ){
      # Check if recombination events and reference genome annotation file is available
      if( !is.null(gubbins.gff.file) & !is.null(ref.genome.name) & is.character(ref.genome.name) ){
        # Check if genome annotation feature labels and genome annotations should be shown
        if( isTRUE(show.gene.label) & isTRUE(show.genome.annot) ){
          # Show genome annotation feature labels
          par(mai=c(0,0,0,0))
          plot(genome.start,1,las=1,xlim=c(genome.start,genome.end ),ylim=c(0,10),bty="n",xaxt="n",yaxt="n",
               xlab="",ylab="",col=rgb(0,0,0,alpha=0),xaxs="i",yaxs="i")
          gb.file2<-load.genome.GFF(ref.genome.name) %>%
            dplyr::filter(.data$start>=as.vector(genome.start) & .data$end<=as.vector(genome.end))
          text((gb.file2$start+gb.file2$end)/2,0.10,
               labels=gb.file2$gene,
               col=ifelse(gb.file2$start>=genome.start & gb.file2$end<=genome.end,"black","black"),
               srt=gene.label.angle,adj=0)
        }else{
          # Hide genome annotation feature labels
          # Empty plot as a placeholder
          par(mai=c(0,0,0,0))
          show.blank.plot()
        }
      }else{
        # Check if genome annotations are provide through a data frame
        if( length(base::setdiff(class(ref.genome.name),c("tbl_df","tbl","grouped_df","data.frame","rowwise_df")))==0  ){
          # Check if genome annotation feature labels should be shown
          if( isTRUE(show.gene.label) ){
            par(mai=c(0,0,0,0))
            plot(genome.start,1,las=1,xlim=c(genome.start,genome.end ),ylim=c(0,10),bty="n",xaxt="n",yaxt="n",
                 xlab="",ylab="",col=rgb(0,0,0,alpha=0),xaxs="i",yaxs="i")
            gb.file2<-ref.genome.name %>% dplyr::filter(.data$feature!="source") %>%
              dplyr::filter(.data$start>=as.vector(genome.start) & .data$end<=as.vector(genome.end))
            text((gb.file2$start+gb.file2$end)/2,0.10,
                 labels=gb.file2$gene,
                 col=ifelse(gb.file2$start>=genome.start & gb.file2$end<=genome.end,"black","black"),
                 srt=gene.label.angle,adj=0)
          }else{
            # Hide genome annotation feature labels
            # Empty plot as a placeholder
            par(mai=c(0,0,0,0))
            show.blank.plot()
          }
        }else{
          # Hide genome annotation feature labels
          # Empty plot as a placeholder
          par(mai=c(0,0,0,0))
          show.blank.plot()
        }
      }
    }else{
      # Hide genome annotation feature labels
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }
  ################################  PANEL 11 - RECOMBINATION HEATMAP/FREQUENCY PLOT  #####################################
  {
    ## Show recombination frequency plot in the diagram
    # Check if the recombination frequency plot should be shown
    if( isTRUE(show.rec.freq.per.base) & isTRUE(show.rec.events) ){
      # Show the recombination frequency plot
      par(mai=c(0,0,0,0))
      if(!isTRUE(rec.events.per.base.as.heatmap)){
        genome.ticks <- pretty(c(genome.start,genome.end))

        temp.vals.fr<-count.rec.events.per.base(gubbins.gff.file=gubbins.gff.file,recom.input.type=recom.input.type)
        plot(0,0,las=1,xlim=c(genome.start,genome.end),ylim=c(0,max(temp.vals.fr$FRQ)+1 ),
             bty="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i",col=rgb(0,0,0,alpha=0),
             xlab="Chromosome position (bp)",ylab=expression("N"[rec]),xaxt="n")
        axis(1,at=genome.ticks,lwd.ticks=show.genome.ticks,lwd=show.genome.axis,
             labels=formatC(pretty(c(genome.start,genome.end)),digits = 0,format = "d" ),las=1)
        axis(2,at=genome.ticks,lwd.ticks=show.genome.ticks,lwd=show.genome.axis, #lwd.ticks=1,lwd=0,
             labels=formatC(pretty(c(genome.start,genome.end)),digits = 0,format = "d" ),las=1)
        polygon(c(min(temp.vals.fr$POS),temp.vals.fr$POS,max(temp.vals.fr$POS)),
                c(min(temp.vals.fr$FRQ),temp.vals.fr$FRQ,min(temp.vals.fr$FRQ)),
                col=rgb(0,0,0,alpha=1),border=FALSE,lwd=20)
      }else{
        temp.vals.fr<-count.rec.events.per.base(gubbins.gff.file=gubbins.gff.file,recom.input.type=recom.input.type)
        plot(0,0,las=1,xlim=c(genome.start,genome.end),ylim=c(0,1.5 ),
             bty="n",xaxt="n",yaxt="n",xaxs="i",yaxs="r",col=rgb(0,0,0,alpha=0),
             xlab="Chromosome position (bp)",ylab=expression("N"[rec]),xaxt="n")

        temp.vals.fr<-temp.vals.fr %>% dplyr::arrange(.data$POS) %>% dplyr::mutate(N=lead(.data$POS)) %>%
          dplyr::mutate(N=ifelse(is.na(.data$N),.data$POS+1,.data$N)) %>%
          dplyr::mutate(RR=.data$N-.data$POS) %>% ungroup() %>%
          dplyr::mutate(Q=data.table::rleid(.data$RR)) %>%
          dplyr::mutate(XX=paste0(.data$Q,".",.data$FRQ)) %>%
          dplyr::mutate(Q1=data.table::rleid(.data$XX)) %>%
          dplyr::mutate(XX1=paste0(.data$Q1,".",.data$XX)) %>% dplyr::group_by(.data$XX1) %>%
          dplyr::mutate(start=dplyr::first(.data$POS),end=dplyr::last(.data$POS)) %>%
          dplyr::select(.data$FRQ,.data$XX1,.data$start,.data$end) %>%
          dplyr::distinct() %>% dplyr::ungroup() %>%
          dplyr::arrange(.data$end) %>% dplyr::select(.data$FRQ,.data$start,.data$end)

        genome.ticks <- pretty(c(genome.start,genome.end))

        axis(1,at=genome.ticks,lwd.ticks=show.genome.ticks,lwd=show.genome.axis, #lwd.ticks=1,lwd=0,
             labels=formatC(pretty(c(genome.start,genome.end)),digits = 0,format = "d" ),las=1)

        rect(genome.start,0.05,
             genome.end,1.0,
             col = viridis::plasma(ceiling(log(max(temp.vals.fr$FRQ),2)) )[1],
             border = NA, lwd = 1 )

        rect(temp.vals.fr$start,0.05,
             temp.vals.fr$end,1.0,
             col = viridis::plasma(ceiling(log(max(temp.vals.fr$FRQ),2)) )[temp.vals.fr$FRQ],
             border = NA, lwd = 1 )
      }
    }else{
      {
        # Hide the recombination frequency plot
        # Empty plot as a placeholder
        par(mai=c(0,0,0,0))
        show.blank.plot()
      }
    }
  }

  ################################  PANEL 12 - FIGURE LEGEND   #####################################
  {
    ## Generate figure legend based on metadata columns
    # Check if metadata column colour strips are plotted
    if( isTRUE(show.metadata.columns) & !is.null(taxon.metadata.file) ){
      # Show legend for each metadata column
      par(mai=c(0,0,0.1,0))
      plot(genome.start,genome.end,las=1,xlim=c(0,length(taxon.metadata.columns)+10 ),ylim=c(0,10),bty="n",xaxt="n",yaxt="n",
           xlab="",ylab="",col=rgb(0,0,0,alpha=0),xaxs="i",yaxs="r")
      loop.val<-1
      loop.val1<-0

      if(is.null(color.tree.tips.by.column) & !is.null(trait.for.ancestral.reconstr)){
        taxon.metadata.columns.id<-c(taxon.metadata.columns.names,
                                     setdiff(trait.for.ancestral.reconstr,taxon.metadata.columns.names))
      }else if(!is.null(color.tree.tips.by.column) & is.null(trait.for.ancestral.reconstr)){
        taxon.metadata.columns.id<-c(taxon.metadata.columns.names,
                                     setdiff(color.tree.tips.by.column,taxon.metadata.columns.names))
      }else if(!is.null(color.tree.tips.by.column) & !is.null(trait.for.ancestral.reconstr)){
        taxon.metadata.columns.id<-c(taxon.metadata.columns.names,
                                     setdiff(c(color.tree.tips.by.column,trait.for.ancestral.reconstr),taxon.metadata.columns.names))
      }else{
        taxon.metadata.columns.id<-taxon.metadata.columns.names
      }

      for(count.val in rev(taxon.metadata.columns.id)){
        if(is.null(taxon.metadata.columns.colors)){
          strips.tmp<-tmp.data.val[order(tmp.data.val$pos),c("pos",count.val)]
          strips.vals<-gsub("^NA$","N/A",sort(unique(unname(unlist(strips.tmp[,count.val])))))
          strips.tmp.cols<-stats::setNames(color.pallette(length(strips.vals)),strips.vals )
          strips.tmp$col<-strips.tmp.cols[ sapply(unname(unlist(strips.tmp[,2])), as.character)  ]
          strips.tmp<-strips.tmp[,-1] %>% unique()
          colnames(strips.tmp)<-c("trait","col")

          if(show.fig.legend){
            legend(0,loop.val*(10/length(taxon.metadata.columns.id)),fill=strips.tmp$col,legend=strips.tmp$trait,
                   cex=0.75,title=count.val,bty="n",bg="transparent",horiz=TRUE,xjust=0,yjust=1)
          }
          loop.val<-loop.val+1
          loop.val1<-loop.val1+1
        }else{
          if(is.null(color.tree.tips.by.column)){
            strips.tmp<-tmp.data.val[order(tmp.data.val$pos),c("pos",count.val)]
            strips.vals<-gsub("^NA$","N/A",sort(unique(unname(unlist(strips.tmp[,count.val])))))
            strips.tmp.cols<-stats::setNames(color.pallette(length(strips.vals)),strips.vals )
            strips.tmp$col<-strips.tmp.cols[ sapply(unname(unlist(strips.tmp[,2])), as.character)  ]
            strips.tmp<-strips.tmp[,-1] %>% unique()
            colnames(strips.tmp)<-c("trait","col")
            loop.val<-loop.val+1
            loop.val1<-loop.val1+1
          }else{
            if(color.tree.tips.by.column==count.val){
              strips.tmp<-tmp.data.val[order(tmp.data.val$pos),c("pos",count.val)]
              strips.vals<-gsub("^NA$","N/A",sort(unique(unname(unlist(strips.tmp[,count.val])))))
              strips.tmp.cols<-stats::setNames(color.pallette(length(strips.vals)),strips.vals )
              strips.tmp$col<-strips.tmp.cols[ sapply(unname(unlist(strips.tmp[,2])), as.character)  ]
              strips.tmp<-strips.tmp[,-1] %>% unique()
              colnames(strips.tmp)<-c("trait","col")
              loop.val1<-loop.val1+1
            }else{
              strips.tmp<-tmp.data.val[order(tmp.data.val$pos),c("pos",count.val)]
              strips.vals<-gsub("^NA$","N/A",sort(unique(unname(unlist(strips.tmp[,count.val])))))
              strips.tmp.cols<-stats::setNames(color.pallette(length(strips.vals)),strips.vals )
              strips.tmp$col<-strips.tmp.cols[ sapply(unname(unlist(strips.tmp[,2])), as.character)  ]
              strips.tmp<-strips.tmp[,-1] %>% unique()
              colnames(strips.tmp)<-c("trait","col")

              loop.val<-loop.val+1
              loop.val1<-loop.val1+1
            }
          }

          if(show.fig.legend){
            legend(0,loop.val1*(10/length(taxon.metadata.columns.id)),fill=strips.tmp$col,legend=strips.tmp$trait,
                   cex=0.75,title=count.val,bty="n",bg="transparent",horiz=TRUE,xjust=0,yjust=1)
          }
        }
      }
    }else{
      # Hide legend for each metadata column
      # Empty plot as a placeholder
      par(mai=c(0,0,0,0))
      show.blank.plot()
    }
  }

  ################################  PANEL 13 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 14 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 15 - RECOMBINATION EVENTS PER GENOME  #####################################
  {
    par(mai=c(0,0,0,0))
    if(isTRUE(show.rec.freq.per.genome) & !is.null(gubbins.gff.file) ){
      if( isTRUE(show.rec.events) ){
        if(isTRUE(show.rec.per.genome.scale)){
          rec.events.per.taxon<-count.rec.events.per.genome(gubbins.gff.file=gubbins.gff.file,taxon.names=names(taxon.names) )
          barplot(rec.events.per.taxon,xaxs="i",yaxs="r",border=NA,
                  horiz=TRUE,yaxt="n",xlab=NA,ylab=NA)
        }else{
          rec.events.per.taxon<-count.rec.events.per.genome(gubbins.gff.file=gubbins.gff.file,taxon.names=names(taxon.names) )
          barplot(rec.events.per.taxon,xaxs="i",yaxs="r",border=NA,
                  horiz=TRUE,yaxt="n",xaxt="n",xlab=NA,ylab=NA)
        }
      }else{
        show.blank.plot()
      }
    }
  }

  ################################  PANEL 16 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 17 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 18 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  ################################  PANEL 19 - PLACEHOLDER  #####################################
  {
    ## Empty plot as a placeholder
    par(mai=c(0,0,0,0))
    show.blank.plot()
  }

  # Check if the plot should be viewed in R or saved to a file
  if( !is.null(save.to.this.file) ){
    invisible(dev.off())
    par(mfrow=c(1,1),mai=c(1,1,1,1))
  }else{
    par(mfrow=c(1,1),mai=c(1,1,1,1))
  }
}

