#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package:  
#

#' Accepts a new gene.
#' 
#' @description Executes a genetic operator pipeline.
#'              The new gene is returned.
#'
#' @param OperatorPipeline   Genetic operator pipeline (an R function). 
#' @param gene               Gene.
#' @param lF                 Local configuration.
#'
#' @return New gene.
#'
#' @family Acceptance Rule
#'
#' @examples
#' id<-function(x, lF){x}
#' g1<-InitGene(lFxegaGaGene)
#' AcceptNewGene(id, g1, lFxegaGaGene)
#' @export
AcceptNewGene<-function(OperatorPipeline, gene, lF)
{  
newGene<-OperatorPipeline(gene, lF)
  if (newGene$evalFail) { 
	  warning("AcceptNewGene: Evaluation failed. Gene: \n", newGene)
  }
return(newGene)}

#' Accepts only genes with equal or better fitness.
#'
#' @description Change the gene by a genetic operator pipeline 
#'              and return the new gene only if the new gene 
#'              has at least the same fitness as the gene.
#'
#' @details  The fitness of genes never decreases.
#'           New genes with inferior fitness do not survive.
#'
#' @param OperatorPipeline    Genetic operator pipeline. 
#' @param gene                Gene.
#' @param lF                  Local configuration.
#'
#' @return The new gene, if it is at least as fit as \code{gene} else
#'         the old gene \code{gene}.
#'
#' @family Acceptance Rule
#'
#' @examples
#' OPpipe1<-function(g, lF){InitGene(lF)}
#' g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
#' g2<-AcceptBest(OPpipe1, g1, lFxegaGaGene)
#' identical(g1, g2)
#' @export
AcceptBest<-function(OperatorPipeline, gene, lF)
{ newGene<- lF$EvalGene(OperatorPipeline(gene, lF), lF)
  if (newGene$evalFail) { 
	  warning("AcceptBest: Evaluation failed. Gene: \n", newGene)
  }
  if(newGene$fit>=gene$fit) 
    {return(newGene)} else {return(gene)}}

#' Metropolis Acceptance Rule. 
#'
#' @description Change the gene by a genetic operator pipeline.
#'              Always accept a new gene with a fitness improvement.
#'              For maximizing fitness
#'              accept genes with lower fitness with probability
#'              \code{(runif(1)<exp(-(fitness-newfitness)*beta/Temperature)}
#'              and reduce temperature with a cooling schedule.
#'              Used: \code{Temperature<-alpha*Temperature} with 
#'              \code{alpha<1}.
#'
#' @details The temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'
#' @param OperatorPipeline    Genetic operator pipeline.
#' @param gene                Gene.
#' @param lF                  Local configuration.
#'
#' @return The new gene if it has at least equal performance as the 
#'         old gene else the old gene.
#'
#' @references 
#'      Kirkpatrick, S., Gelatt, C. D. J, and Vecchi, M. P. (1983):
#'      Optimization by Simulated Annealing. 
#'      Science, 220(4598): 671-680.
#'      <doi:10.1126/science.220.4598.671>
#'
#'      Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H.,
#'      Teller, E. (1953):
#'      Equation of state calculations by fast computing machines.
#'      Journal of Chemical Physics, 21(6):1087 – 1092.
#'      <doi:10.1063/1.1699114> 
#'     
#' @family Acceptance Rule
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lFxegaGaGene$Beta<-parm(1)
#' lFxegaGaGene$TempK<-parm(10)
#' OPpipe1<-function(g, lF){InitGene(lF)}
#' g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
#' g2<-AcceptMetropolis(OPpipe1, g1, lFxegaGaGene)

#' @importFrom stats runif    
#' @export
AcceptMetropolis<-function(OperatorPipeline, gene, lF)
{ newGene<- lF$EvalGene(OperatorPipeline(gene, lF), lF)
  if (newGene$evalFail) { 
	  warning("Accept Best: Evaluation failed. Gene: \n", newGene)
	 }
  if(newGene$fit>=gene$fit) {return(newGene)}
  d<-gene$fit-newGene$fit
  if (runif(1) < exp(-d*lF$Beta()/lF$TempK()))  
  {return(newGene)}
  else {return(gene)}
   }

#' Individually Adaptive Metropolis Acceptance Rule. 
#'
#' @description Change the gene by a genetic operator pipeline.
#'              Always accept new genes with fitness improvement.
#'              For maximizing fitness
#'              accept genes with lower fitness with probability
#'              \code{(runif(1)<exp(-(fitness-newfitness)*beta/Temperature)}
#'              and reduce temperature with a cooling schedule.
#'              For each gene, the temperature is corrected upward by a term 
#'              whose size is proportional to 
#'              the difference between the fitness of the current best gene 
#'              in the population and the fitness of the gene.
#'
#' @details The temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'
#' @param OperatorPipeline    Genetic operator pipeline.
#' @param gene                Gene.
#' @param lF                  Local configuration. 
#'
#' @return The new gene if it has at least equal performance as the old gene 
#'         else the old gene.
#'
#' @references 
#'      Locatelli, M. (2000):
#'      Convergence of a Simulated Annealing Algorithm for 
#'      Continuous Global Optimization.
#'      Journal of Global Optimization, 18:219-233.
#'      <doi:10.1023/A:1008339019740>
#'
#' The-Crankshaft Publishing (2023):
#' A Comparison of Cooling Schedules for Simulated Annealing.
#' <URL:https://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/> 
#'
#' @family Acceptance Rule
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lFxegaGaGene$Beta<-parm(1)
#' lFxegaGaGene$TempK<-parm(10)
#' set.seed(2)
#' OPpipe1<-function(g, lF){InitGene(lF)}
#' g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
#' lFxegaGaGene$CBestFitness<-parm(g1$fit)
#' g2<-AcceptMetropolis(OPpipe1, g1, lFxegaGaGene)
#' @importFrom stats runif    
#' @export
AcceptIVMetropolis<-function(OperatorPipeline, gene, lF)
{ newGene<- lF$EvalGene(OperatorPipeline(gene, lF), lF)
  if (newGene$evalFail) { 
	  warning("Accept Best: Evaluation failed. Gene: \n", newGene)
	  }
  if(newGene$fit>=gene$fit) {return(newGene)}
  d<-gene$fit-newGene$fit
  # Individually Variable adaptive correction of temperature:
  Tcorrected<-lF$TempK()*min(2,(1+((lF$CBestFitness()-gene$fit)/gene$fit)))
  if (runif(1) < exp(-d*lF$Beta()/Tcorrected))  
  {return(newGene)}
  else {return(gene)}
   }

#' Metropolis acceptance probability.
#'
#' @param d            Distance between the fitness of the old and the new gene.
#' @param beta         Constant.
#' @param temperature  Temperature.
#'
#' @return  Acceptance probability.  
#'
#' @family Diagnostic
#'
#' @examples
#' MetropolisAcceptanceProbability(d=0, beta=1, temperature=10)
#' MetropolisAcceptanceProbability(d=1, beta=1, temperature=10)
#'@export
MetropolisAcceptanceProbability<-function(d, beta, temperature)
{
	exp(-d*beta/temperature)
}

#' Metropolis acceptance probability table.
#'
#' @param d            Distance between the fitness of the old and the new gene.
#' @param beta         Constant.
#' @param temperature  Temperature.
#' @param alpha        Cooling constant in [0, 1].  
#' @param steps        Number of steps.
#'
#' @return Data frame with the columns alpha, beta, temperature, 
#'         d (distance between fitness), and probability of acceptance.
#'
#' @family Diagnostic
#'
#' @examples
#' MetropolisTable(d=2, beta=2, temperature=10, alpha=0.99, steps=10)
#'@export
MetropolisTable<-function(d=1, beta=2, temperature=1000, alpha=0.9, steps=1000)
{
	df<-data.frame(matrix(ncol=5, nrow=steps))
	colnames(df)<-c("alpha", "beta", "Temperature", "d", "P(Accept)")
        for (i in (1:steps))
	{
          df[i,]<-c(alpha, beta, temperature, d, 
		    MetropolisAcceptanceProbability(d, beta, temperature))
              temperature<-alpha*temperature	
	}
	return(df)
}


#' Configure the acceptance function of a genetic algorithm.
#'
#' @description \code{AcceptanceFactory()} implements selection
#'              of an acceptance rule.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "All" returns \code{AcceptNewGene()} (Default).
#'              \item "Best" returns \code{AcceptBest()}.
#'              \item "Metropolis" returns \code{AcceptMetropolis()}.
#'              \item "IVMetropolis" returns \code{AcceptIVMetropolis()}.
#'              }
#'
#' @param method   A string specifying the acceptance rule.
#'
#' @return An acceptance rule for genes.
#'
#' @family Configuration
#'
#' @export
AcceptFactory<-function(method="All") {
if (method=="All") {f<- AcceptNewGene}
if (method=="Best") {f<- AcceptBest}
if (method=="Metropolis") {f<- AcceptMetropolis}
if (method=="IVMetropolis") {f<- AcceptIVMetropolis}
if (!exists("f", inherits=FALSE))
        {stop("Acceptance label ", method, " does not exist")}
return(f)
}

