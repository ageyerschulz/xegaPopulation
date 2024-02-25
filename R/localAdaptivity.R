#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package:  
#

#' Individually adaptive mutation rate. (Bit mutation Rate)
#'
#' @description Adaptivity of a local operator mutation parameter.
#'              Currently not used. Implements a threshold rule.
#'              The rule is implemented directly in IVAdaptiveMutateGene.
#'              in package xegaGaGene. Move? 
#'
#' @details     TODO: Move this xegaGaGene and generalize the
#'              bit mutation operator and introduce a factory 
#'              for bit mutation rates. Rationale: Local parameters
#'              are representation dependent.      
#'
#' @param fit  Fitness of gene.
#' @param lF   Local configuration.
#'
#' @return Mutation rate of a gene depending on its fitness.
#'  
#' @export
IAMBitRate<-function(fit, lF)
{
if (fit>(lF$CutoffFit()*lF$CBestFitness()))
                {lF$BitMutationRate1()}
        else
                {lF$BitMutationRate2()}
}


