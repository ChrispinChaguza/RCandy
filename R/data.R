# This goes in R/data.R

#' @title Gubbins recombination data
#' @description This dataset provides an example phylogenetic tree, strain metadata, reference genome and recombination events for 170 Streptococcus pneumoniae isolates collected globally belonging to the multilocus sequence typing (MLST) clone ST320.
#' @format A list containing four data frames:
#' \describe{
#'   \item{\code{tree}}{phylo: Phylogenetic tree in ape's phylo format}
#'   \item{\code{metadata}}{data frame: Data frame containing metadata for the strains in the tree}
#'   \item{\code{gubbins.GFF}}{data frame: Data frame containing recombination events identified in each strain (GFF format)}
#'   \item{\code{refgenome.GFF}}{data frame: Data frame containing genomic features in the reference genome (GFF format)}
#'}
#' @docType data
#' @keywords datasets
#' @name "RCandy"
#' @usage data("RCandy")
#' @format A list containing pre-loaded phylogenetic tree, taxon metadata, Gubbins GFF file, and reference genome GFF file
#'
#' @source Gladstone RA, Lo SW et al. International genomic definition of pneumococcal lineages, to contextualise disease, antibiotic resistance and vaccine impact. EBioMedicine. 2019 May;43:338-346. doi: 10.1016/j.ebiom.2019.04.021. Epub 2019 Apr 16. PMID: 31003929; PMCID: PMC6557916.
#'
#' @source The Global Pneumococcal Sequencing (GPS) Consortium: https://www.pneumogen.net/gps/
#'
#' @source The reference whole genome sequence was obtained from GenBank (accession number: NC_010380): https://www.ncbi.nlm.nih.gov/nuccore/NC_010380.1.

"RCandy"
