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
#' @details The current version is sequential.
#'          For parallelization, a restructuring of the 
#'          main loop with an integration of \code{xegaNextPopulation}
#'          and \code{xegaEvalPopulation} is planned, because
#'          this allows the parallelization of a large part of 
#'          the genetic operations which are sequential in the 
#'          current version. 
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
newlF$CBestFitness<-xegaSelectGene::parm(max(fit))
newlF$CMeanFitness<-xegaSelectGene::parm(mean(fit))
newlF$CVarFitness<-xegaSelectGene::parm(var(fit))
newlF$CWorstFitness<-xegaSelectGene::parm(min(fit))

if (lF$SelectionContinuation()==TRUE) 
{ newlF$SelectGene<-TransformSelect(fit, lF, lF$SelectGene)
 newlF$SelectMate<-TransformSelect(fit, lF, lF$SelectMate)}

while(length(newpop)<length(pop))
 {
newpop<-append(newpop, lF$ReplicateGene(pop, fit, newlF))
 }
return(head(newpop,length(pop)))
}

