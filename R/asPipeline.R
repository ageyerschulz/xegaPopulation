
# as Pipeline 
# (c) 2025 Andreas Geyer-Schulz.

#' Identity (No compilation of genetic operator pipelines for population).
#'
#' @param pop  A population.
#' @param lF   The local function configuration.
#'
#' @return The population
#'
#' @family Genetic operator pipelines
#'
#' @examples
#' pop5<-xegaInitPopulation(5, lFxegaGaGene)
#' pop5c<-asPipeline(pop5, lFxegaGaGene)
#' identical(pop5, pop5c)
#' @export
asPipelineID<-function(pop, lF) {pop}

#' Converts a population into a list of genetic operator pipelines.
#'
#' @param pop  A population.
#' @param lF   The local function configuration.
#'
#' @return A list of genetic operator pipelines (closures).
#'
#' @family Genetic operator pipelines
#'
#' @examples
#' pop5<-xegaInitPopulation(5, lFxegaGaGene)
#' pop5c<-asPipeline(pop5, lFxegaGaGene)
#' @importFrom xegaGaGene newPipeline
#' @export
asPipeline<-function(pop, lF)
{unlist(lapply(pop, xegaGaGene::newPipeline, lF=lF))}

#' Embeds genetic operator pipelines into the genes of a population.
#'
#' @param pop  A population.
#' @param lF   The local function configuration.
#'
#' @return A population of genes with embedded genetic operator pipelines.
#'
#' @family Genetic operator pipelines
#'
#' @examples
#' pop5<-xegaInitPopulation(5, lFxegaGaGene)
#' pop5c<-asPipelineG(pop5, lFxegaGaGene)
#' @importFrom xegaGaGene newPipelineG
#' @export
asPipelineG<-function(pop, lF)
{lapply(pop, xegaGaGene::newPipelineG)}

#' Configure asPipeline. 
#'
#' @description \code{xegaAsPipelineFactory()} implements the selection
#'              of one of the asPipeline functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "NoPipe" returns \code{asPipelineID()}:
#'                    No change to population.
#'              \item "PipeC" returns \code{asPipeline()}: 
#'                    Population a list of closures of genetic operator 
#'                    populations.
#'              \item "PipeG" returns \code{asPipelineG()}:
#'                    Population consists of genes with embedded 
#'                    genetic operator pipelines.
#'              }
#'
#' @param method     A string specifying the asPipeline function.
#'
#' @return An asPipeline function.
#'
#' @family Configuration
#'
#' @examples
#' xegaAsPipelineFactory("PipeC")
#' @export
xegaAsPipelineFactory<-function(method="NoPipe") {
if (method=="NoPipe") {f<- asPipelineID}
if (method=="PipeC") {f<- asPipeline}
if (method=="PipeG") {f<- asPipelineG}
if (!exists("f", inherits=FALSE))
        {stop("population asPipeline label ", method, " does not exist")}
return(f)
}

