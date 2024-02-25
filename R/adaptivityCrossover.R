#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package:  
#

#' Individually adaptive crossover rate.
#'
#' @description The basic idea is to apply crossover to a gene whose
#'              fitness is below a threshold value with higher probability 
#'              to give it a chance 
#'              to improve. The threshold value is computed by
#'              \code{lF$CutoffFit()*lF$CBestFitness()}. 
#'
#' @details The following constants are used:
#'         \code{lF$CrossRate1()<lF$CrossRate2()}, and 
#'         \code{lF$CutoffFit()} in [0, 1].
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
#' @return Crossover rate of a gene depending on its fitness.
#'  
#' @family Adaptive Rates
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$CrossRate1<-parm(0.20) 
#' lF$CrossRate2<-parm(0.40) 
#' lF$CutoffFit<-parm(0.60)
#' lF$CBestFitness<-parm(105)
#' IACRate(100, lF)
#' IACRate(50, lF)
#' @export
IACRate<-function(fit, lF)
{
if (fit>(lF$CutoffFit()*lF$CBestFitness()))
                {lF$CrossRate1()}
        else
                {lF$CrossRate2()}
}

#' Constant crossover rate.
#'        
#'
#' @param fit  Fitness of gene.
#' @param lF   Local configuration.
#'
#' @return Constant crossover rate.
#'  
#' @family Rates
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(CrossRate1=parm(0.20))
#' ConstCRate(100, lF)
#' ConstCRate(50, lF)
#' @export
ConstCRate<-function(fit, lF)
{ lF$CrossRate1() }

#' Configure the crossover function of a genetic algorithm.
#'
#' @description \code{CrossRateFactory()} implements selection
#'              of one of the crossover rate functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Const" returns \code{ConstCRate()}.
#'              \item "IV" returns \code{IACrate()}.
#'                    This function gives bad genes a higher cross rate.
#'              }
#'
#' @param method A string specifying a function for the crossover rate.
#'
#' @return Crossover rate function.
#'
#' @family Configuration
#'
#' @examples
#' f<-CrossRateFactory("Const")
#' f(10, list(CrossRate1=function() {0.2}))
#' @export
CrossRateFactory<-function(method="Const") {
if (method=="Const") {f<- ConstCRate}
if (method=="IV") {f<- IACRate}
if (!exists("f", inherits=FALSE))
        {stop("Crossover Rate label ", method, " does not exist")}
return(f)
}

