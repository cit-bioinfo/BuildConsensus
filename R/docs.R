#' Subtype annotation for 1750 BLCA samples in 6 different classification systems
#' 
#' A dataset containing class annotation for muscle-invasive bladder 
#' cancer samples in 6 different published classification systems.
#' 
#' @format A data frame with 1750 observations of 6 variables:
#' \describe{
#'    \item{Baylor.subtype}{A subtyping system consisting of two classes: 
#'    'Basal' and 'Differentiated'}
#'    \item{ChapelHill.subtype}{A subtyping system consisting of two classes: 
#'     'Luminal' and 'Basal'.}
#'     \item{CIT.subtype}{A subtyping system consisting of six classes:
#'     MC1 to MC6}
#'     \item{MDA.subtype}{A subtyping system consisting of three classes: 
#'    'luminal', 'p53-like' and 'basal'}
#'    \item{TCGA2017.subtype}{A subtyping system consisting of five classes: 
#'    'Luminal_papillary', 'Luminal_infiltrated', 'Luminal', Neuronal' and 'Basal_squamous'}
#'    \item{Lund2017.subtype}{A subtyping system consisting of ten classes: 
#'    'GU', 'GU-Inf', 'UroA-Prog', 'Uro-Inf', 'UroB', 'UroC', 'Mes-like', 'Ba/Sq', 
#'    'Ba/Sq-Inf','Sc/NE-like'}
#' }
#' 
#' @references \url{https://www.biorxiv.org/content/10.1101/488460v2}
#' @author Aurelie Kamoun
"blca_class"

#' Build a Consensus Classification from Multiple Classification Systems
#' 
#' This package implements methods for building a consensus molecular 
#' classification using existing classification systems.
#' 
#' @seealso \link[BuildConsensus]{consensus.MCL}, \link[BuildConsensus]{consensus.CIT},
#' \link[BuildConsensus]{consensus.COCA}
#' @author Aurelien de Reynies, Aurelie Kamoun
#' @note This is a contribution from the Tumor Identity Cards (CIT) program founded by the 
#' 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. 
#' For any question please contact \url{CITR@ligue-cancer.net}
#' 
"_PACKAGE"