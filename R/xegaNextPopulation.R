#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Population-level functions.
#                 Independent of gene representation.
#                 The replication mechanism and its variants
#          Package: xegaPopulation.
#

#' Import for examples.
#' @importFrom xegaGaGene lFxegaGaGene
#' @export
lFxegaGaGene<-xegaGaGene::lFxegaGaGene

#' Import for examples.
#' @param lF   a list of local functions
#' @return a new random gene
#' @importFrom xegaGaGene xegaGaInitGene
#' @export
InitGene<-xegaGaGene::xegaGaInitGene

#' Import for examples.
#'
#' @param gg1 a gene
#' @param gg2 a gene
#' @param lF  list of local functions
#'
#' @return a list of one gene
#' 
#' @importFrom xegaGaGene  xegaGaCrossGene
#' @export 
CrossGene<-xegaGaGene::xegaGaCrossGene

#' Import for examples.
#'
#' @param gg1 a gene
#' @param gg2 a gene
#' @param lF  list of local functions
#'
#' @return a list of two genes
#' 
#' @importFrom xegaGaGene  xegaGaCross2Gene
#' @export 
Cross2Gene<-xegaGaGene::xegaGaCross2Gene

#' Import for examples.
#'
#' @param pop the population.
#' @param fit the fitness-
#' @param lF  list of local functions
#'
#' @return a list with one gene
#' 
#' @importFrom xegaGaGene xegaGaReplicateGene
#' @export 
ReplicateGene<-xegaGaGene::xegaGaReplicateGene


#' Computes the next population of genes.
#'
#' @description \code{xegaNextPopulation()} 
#'              builds the next population by repeatedly
#'              calling \code{ReplicateGene()}.
#'
#' @details Generating the next population is sequential. 
#'          However, in order to shift more computations 
#'          into the evaluation step, genetic operator pipelines
#'          have been implemented by 
#'          \code{lF$ReplicateGene()}. 
#'          \code{xegaNextPopulation()} only has to convert 
#'          elitist solutions into function closures.            
#'
#'          For adaptive genetic operators, population statistics and 
#'          the current generation are stored as constant functions 
#'          in \code{lF}.
#'
#' @param pop     Population of genes.
#' @param fit     Fitness.
#' @param lF      Local configuration.
#'
#' @return Population of genes.
#'
#' @family Population Layer
#'
#' @examples
#' lFxegaGaGene$cGeneration<-function() {0}
#' lFxegaGaGene$MutationRate<-MutationRateFactory(method="Const")
#' lFxegaGaGene$ReplicateGene<-ReplicateGene
#' lFxegaGaGene$Accept<-AcceptFactory(method="All")
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' newpop<-xegaNextPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
#'
#' @importFrom stats         var
#' @importFrom xegaSelectGene  TransformSelect
#' @importFrom xegaSelectGene  parm
#' @export
xegaNextPopulation<-function(pop, fit, lF)
{ 
if (lF$Elitist()) {
   newpop<-list(pop[[xegaBestGeneInPopulation(fit)[1]]])
   if (lF$Pipeline()==TRUE) 
     {newpop<-asPipeline(newpop, lF)}
   } else {
   newpop<-list() }

newlF<-lF

# cat("xegaNextPopulation:\n")
# cat("CBestFitness:", newlF$CBestFitness(), "\n")
# cat("cGeneration:", newlF$cGeneration(), "\n")

if (lF$SelectionContinuation()==TRUE) 
{ newlF$SelectGene<-TransformSelect(fit, lF, lF$SelectGene)
 newlF$SelectMate<-TransformSelect(fit, lF, lF$SelectMate)}

while(length(newpop)<length(pop))
 {
newpop<-append(newpop, lF$ReplicateGene(pop, fit, newlF))
 }
return(head(newpop,length(pop)))
}

