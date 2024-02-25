#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package:  
#

#' Individually adaptive mutation rate.
#'
#' @description The probability of applying a mutation operator
#'              to a gene. The idea is that a gene selected for 
#'              reproduction whose fitness is 
#'              below a threshold value is mutated with a higher 
#'              probability to give it a chance.  
#'
#' @details The probability of applying a mutation operator is
#'          determined by a threshold: If the fitness of a gene
#'          is higher than \code{lF$CutoffFit()*lF$CBestFitness()},
#'          than return \code{lF$MutationRate1()} 
#'          else \code{lF$MutationRate2()}.
#'
#'          Note that the idea is also applicable to gene specific 
#'          local mutation operators. For example, the bit mutation rate
#'          of mutation operators for binary genes.
#'
#' @references
#'   Stanhope, Stephen A. and Daida, Jason M. (1996)
#'   An Individually Variable Mutation-rate Strategy for Genetic Algorithms.
#'   In: Koza, John (Ed.)
#'   Late Breaking Papers at the Genetic Programming 1996 Conference.
#'   Stanford University Bookstore, Stanford, pp. 177-185.
#'   (ISBN:0-18-201-031-7)
#'
#' @param fit  Fitness of gene.
#' @param lF   Local configuration.
#'
#' @return Mutation rate of a gene depending on its fitness.
#'  
#' @family Adaptive Rates
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$MutationRate1<-parm(0.20)
#' lF$MutationRate2<-parm(0.40)
#' lF$CutoffFit<-parm(0.60)
#' lF$CBestFitness=parm(105)
#' IAMRate(100, lF)
#' IAMRate(50, lF)
#' @export
IAMRate<-function(fit, lF)
{
if (fit>(lF$CutoffFit()*lF$CBestFitness()))
                {lF$MutationRate1()}
        else
                {lF$MutationRate2()}
}

#' Constant mutation rate.
#'        
#' @param fit  Fitness of gene.
#' @param lF   Local configuration.
#'
#' @return Constant mutation rate.
#'  
#' @family Rates
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$MutationRate1<-parm(0.20)
#' ConstMRate(100, lF)
#' ConstMRate(50, lF)
#' @export
ConstMRate<-function(fit, lF)
{ lF$MutationRate1() }

#' Configure the mutation rate function of a genetic algorithm.
#'
#' @description The \code{MutationRateFactory()} implements selection
#'              of one of the crossover rate functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Const" returns \code{ConstMRate()} (Default).
#'              \item "IV" returns \code{IAMrate()}.
#'                    This function gives bad genes a higher mutation rate.
#'              }
#'
#' @param method   A string specifying a function for the mutation rate.
#'
#' @return A mutation rate function.
#'
#' @family Configuration
#'
#' @examples
#' f<-MutationRateFactory("Const")
#' f(10, list(MutationRate1=function() {0.2}))
#' @export
MutationRateFactory<-function(method="Const") {
if (method=="Const") {f<- ConstMRate}
if (method=="IV") {f<- IAMRate}
if (!exists("f", inherits=FALSE))
        {stop("Mutation Rate label ", method, " does not exist")}
return(f)
}

