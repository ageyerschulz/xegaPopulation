#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Population-level functions.
#                 Independent of gene representation.
#          Package: xegaPopulation.
#

#' Import lFxegaGaGene
#' @importFrom xegaGaGene lFxegaGaGene
lFxegaGaGene<-xegaGaGene::lFxegaGaGene

#' Initializes a population of genes.
#'
#' @description \code{xegaInitPopulation()} initializes a population
#'                 of genes.
#'
#' @param popsize    Population size.
#' @param lF         Local function configuration.
#' 
#' @return List of genes.
#'
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#'
#' @export
xegaInitPopulation<-function(popsize, lF)
{ pop<-list()
  for (i in 1:popsize)
  { pop[[i]]<-lF$InitGene(lF)}
  return(pop) }

#' Observe summary statistics of the fitness of the population.
#'
#' @description \code{xegaObservePopulation()} reports
#'              summary statistics of the fitness of the population.
#'
#' @details Population statistics are used for 
#'              \itemize{
#'                 \item implementing individually variable operator rates and
#'                 \item visualizing the progress of the algorithm.      
#'                 }
#'
#' @param fit     Vector of fitness values of a population.
#' @param v       Vector of population statistic vectors.
#'
#' @return Vector of population statistics. If position
#'         \code{x} modulo \code{8} equals
#'         \enumerate{
#'         \item \code{1}:  Mean fitness.
#'         \item \code{2}:  Min fitness.
#'         \item \code{3}:  Lower-hinge 
#'                        (approx. 1st quartile) of fitness.
#'         \item \code{4}:  Median fitness.
#'         \item \code{5}:  Upper-hinge 
#'                         (approx. 3rd quartile) of fitness.
#'         \item \code{6}:  Max fitness.
#'         \item \code{7}:  Variance.
#'         \item \code{8}: Mean absolute deviation. 
#'             }
#'
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' popStats<-xegaObservePopulation(epop10$fit)
#' popStats<-xegaObservePopulation(epop10$fit, popStats)
#' matrix(popStats, ncol=8, byrow=TRUE)
#'
#' @importFrom stats fivenum
#' @importFrom stats mad
#' @export
xegaObservePopulation<-function(fit, v=vector())
{
return(append(v, c(mean(fit), fivenum(fit), var(fit), mad(fit, constant=1))))
}

#' Combine fitness, generations, and the phenotype of the gene. 
#' 
#' @param pop          Population.
#' @param evallog      Evaluation log.
#' @param generation   Generation logged.
#' @param lF           Local function configuration.
#'
#' @return Update of the evaluation log. 
#'         The evaluation log is a list of decoded and evaluated genes.
#'         A list item of the evaluation log has the following 
#'         elements:
#'         \itemize{
#'         \item \code{$generation}:   The generation.
#'         \item \code{$fit}:          The fitness value.
#'         \item \code{$sigma}:        The standard deviation of the 
#'                                     fitness value, if it exists.
#'                                     Default: \code{0}.
#'         \item \code{$obs}:          The number of observations for 
#'                                     computing the 
#'                                     fitness value, if it exists.
#'                                     Default: \code{0}.
#'         \item \code{$phenotype}:    The phenotype of the gene.
#'         }
#'
#' @family Population Layer
#'
#'@examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' logevals<-list()
#' logevals
#' logevals<-xegaLogEvalsPopulation(epop10$pop, logevals, 1, lFxegaGaGene)
#' logevals
#'@export
xegaLogEvalsPopulation<-function(pop, evallog, generation, lF)
{
	new<-list()
        sigmaExists<-("sigma" %in% names(pop$pop[[1]]$sigma))
	for (i in (1:length(pop))) {
        v<-list()
        v$generation<-generation
        v$fit<-pop[[i]]$fit
        v$sigma<-v$obs<-0
        if (sigmaExists)  
            {v$sigma<-pop[[i]]$sigma
             v$obs<-pop[[i]]$obs}
        v$phenotype<-lF$DecodeGene(pop[[i]], lF)
	new[[i]]<-v}
	return(c(evallog,new))
}

#' Extracts indices of best genes in population.
#'
#' @description \code{xegaBestGeneInPopulation()} extracts the indices of
#'              the best genes in the population.
#'
#' @details You might use:   
#'	 \code{which(max(fit)==fit)}. But this is slower! 
#'
#' @param fit   Fitness vector of a population of genes.
#'
#' @return List of the indices of the best genes in the population.
#'
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' xegaBestGeneInPopulation(epop10$fit)
#'
#' @export
xegaBestGeneInPopulation<-function(fit){
        (1:length(fit))[max(fit)==fit]
}

#' Best solution in the population.
#'
#' @description \code{xegaBestInPopulation()} extracts the best
#'              individual of a population and
#'              reports fitness, value, genotype, and phenotype:
#'
#'              \enumerate{
#'              \item
#'              \code{fitness}:  The fitness value of the genetic algorithm.
#'              \item
#'              \code{value}:   The function value of the problem environment.
#'              \item
#'              \code{genotype}:  The gene representation.
#'              \item  
#'              \code{phenotype}: The problem representation.
#'                 E.g. a parameter list, a program, ...
#'              }
#'
#'              We report one of the best solutions.
#'
#' @param pop     Population of genes.
#' @param fit     Vector of fitness values of \code{pop}.
#' @param lF      Local function configuration.
#' @param allsolutions  If TRUE, also return a list of all solutions.
#'
#' @return Named list with the following elements:
#'     \itemize{
#'         \item \code{$name}:     The name of the problem environment.
#'         \item \code{$fitness}:  The fitness value of the best solution.
#'         \item \code{$value}:      The evaluated best gene. 
#'         \item \code{$numberOfSolutions}:   The number of solutions.
#'         \item \code{$genotype}:    The best gene.
#'         \item \code{$phenotype}:   The parameters of the solution
#'                            (the decoded gene).
#'         \item \code{$phenotypeValue}:  The value of the
#'                            function of the parameters of the solution
#'                            (the decoded gene).
#'         \item \code{$allgenotypes}:  The genotypes of all best solutions.
#'                                   (allsolutions==TRUE)
#'         \item \code{$allphenotypes}:  The phenotypes of all best solutions.
#'                                   (allsolutions==TRUE)
#'     }
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' xegaBestInPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
#'
#' @importFrom utils head
#' @export
xegaBestInPopulation<-function(pop, fit, lF, allsolutions=FALSE)
{
        best<-(1:length(fit))[max(fit)==fit]
        bestGene<- pop[[head(best,1)]]
        parms<- lF$DecodeGene(bestGene,lF)
        solution<-list(
             name=list(lF$penv$name()),		       
             fitness=bestGene$fit,
             value=bestGene,
	     numberOfSolutions=length(best),
             genotype=bestGene,
             phenotype=parms, 
	     phenotypeValue=lF$Max()*bestGene$fit)
	if ((allsolutions) && (length(best)>1)) 
	{
	solution[["allgenotypes"]]<-pop[best]
        solution[["allphenotypes"]]<-
		lapply(pop[best], lF$DecodeGene, lF=lF) 	
	}
        return(solution)
}

#' Provide elementary summary statistics of the fitness of the population.
#'
#' @description \code{SummaryPopulation()} reports
#'              on the fitness and the value of the best solution
#'              in the population.
#'              
#'              The value of \code{lF$Verbose()} controls the
#'              information displayed:
#'              \enumerate{ 
#'              \item \code{== 0}: Nothing is displayed.
#'
#'              \item \code{== 1}: 1 point per generation.
#'
#'              \item \code{> 1}: Max(fit), number of solutions, indices.
#'
#'              \item \code{> 2}: and population fitness statistics.
#'       
#'              \item \code{> 3}: and 1 solution.
#'              }
#'
#' @param pop    Population of genes.
#' @param fit    Vector of fitness values of \code{pop}.
#' @param lF     Local function configuration.
#' @param iter   The generation. Default: \code{0}.
#'
#' @return The number \code{0}.
#'
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
#' rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, iter=12)
#'
#' @importFrom utils str
#' @export
xegaSummaryPopulation<-function(pop, fit, lF, iter=0)
{
	if (lF$Verbose()==0) {return(0)}
	if (lF$Verbose()==1) {
		if (0==(iter%%50)) {cat("\n")}
		cat(".")
		return(0)}

	if (lF$Verbose()>1) {
        if (iter==0)
        {
        cat("Best solution:\n")
        } else
        {
        cat("Best Solution in Iteration:", iter, "\n")
        }

        best<-(1:length(fit))[max(fit)==fit]
        cat("Max(fit): ", max(fit),
            "No. solutions: ", length(best),
            "Indices of best genes: ", best,
            "\n")
	}

	if (lF$Verbose()>2) {
	stats<-xegaObservePopulation(fit)
	cat(
	    "Fitness Min:", stats[2],
	    " Q1:", stats[3],
	    " Mean:", stats[1],
	    " Q3:", stats[5],
	    " Max:", stats[6],
            "\n")
	}
	    
	if (lF$Verbose()>3) {
	#cat("In Summary:Best in Population.\n")
        solution<-xegaBestInPopulation(pop, fit, lF)
	cat("Fitness:", solution$fitness, 
	    "Value of Phenotype:", 
	    solution$phenotypeValue, "\n")
	cat("Phenotype:\n")
	print(solution$phenotype)
	}

	if (lF$Verbose()>4) {
	#cat("In Summary:Best in Population.\n")
        solution<-xegaBestInPopulation(pop, fit, lF)
	cat("str(Genotype):\n")
	str(solution$genotype)
	}
        return(0)
}

