
#' Evaluates a population of genes in a problem environment
#'
#' @description \code{xegaEvalPopulation()} evaluates a population
#'                  of genes in a problem environment.
#'
#' @details Parallelization of the evaluation of fitness functions
#'          is possible by defining \code{lF$lapply}.
#'
#' @param pop    Population of genes.
#' @param lF     Local function configuration.
#'
#' @return List of
#'         \itemize{
#'         \item \code{$pop} gene vector,
#'         \item \code{$fit} fitness vector, 
#'         \item \code{$evalFail} number of failed evaluations.
#'         } 
#'
#' @family Population Layer
#'
#' @examples
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' lFxegaGaGene[["lapply"]]<-ApplyFactory(method="Sequential") 
#' result<-xegaEvalPopulation(pop10, lFxegaGaGene)
#'
#' @export
xegaEvalPopulation<-function(pop, lF)
{ pop<- lF$lapply(pop, lF$EvalGene, lF=lF)
  fit<- unlist(lapply(pop, function(x) {x$fit}))
  evalFail<-sum(unlist(lapply(pop, function(x) {x$evalFail})))
  return(list(pop=pop, fit=fit, evalFail=evalFail))
}

#' Evaluates a population of genes in a a problem environment repeatedly.
#'
#' @description \code{xegaRepEvalPopulation()} evaluates a population
#'              of genes in a problem environment \code{lF$rep} times.
#'              The results of repeatedly evaluating a gene are aggregated:
#'              \itemize{
#'              \item \code{gene$fit} is the mean fitness, 
#'              \item \code{gene$var} is the fitness variance, 
#'              \item \code{gene$std} is the standard deviation of the fitness,
#'                    and 
#'              \item \code{gene$obs} is the number of repetitions. 
#'              }
#'
#' @details Parallelization of the evaluation of fitness functions
#'          is possible by defining \code{lF$lapply}.
#'
#' @param pop    Population of genes.
#' @param lF     Local function configuration.
#'
#' @return List of
#'         \itemize{
#'         \item \code{$pop} gene vector,
#'         \item \code{$fit} fitness vector, 
#'         \item \code{$evalFail} number of failed evaluations.
#'         } 
#'
#' @family Population Layer
#'
#' @examples
#'     parm<-function(x){function() {return(x)}}
#' pop10<-xegaInitPopulation(10, lFxegaGaGene)
#' lFxegaGaGene[["lapply"]]<-ApplyFactory(method="Sequential") 
#' lFxegaGaGene[["rep"]]<-parm(3) 
#' result<-xegaRepEvalPopulation(pop10, lFxegaGaGene)
#'
#' @importFrom stats sd
#' @export
xegaRepEvalPopulation<-function(pop, lF)
{ newpop<-unlist(lapply((1:length(pop)), FUN=function(x, p, lF)
                 {rep(p[x], lF$rep())}, p=pop,lF=lF), recursive=FALSE)
index<-lapply((0:(length(pop)-1))*lF$rep(),
       FUN=function(x, lF) {y<-0:(lF$rep()-1); 1+x+y}, lF)
rndPerm<-sample(length(newpop), length(newpop)) # random permutation
rndpop<-newpop[rndPerm]                         # permute
retpop<-lF$lapply(rndpop, lF$EvalGene, lF=lF)      # here: parallel.
retpop[rndPerm]<-retpop                         # inverse permute
tmppop<-pop                                     # begin aggregation
for (i in (1:length(pop)))
    { newtmp<-retpop[index[[i]]]
      new<-newtmp[[1]]
      fit<-(unlist(lapply(newtmp, FUN=function(x) {x$fit})))
      new$fit<-mean(fit); new$var<-var(fit); new$sigma<-sd(fit) 
      new$obs<-length(fit)
      tmppop[[i]]<-new}                         # end aggregation
fit<- unlist(lapply(tmppop, function(x) {x$fit}))
evalFail<-sum(unlist(lapply(tmppop, function(x) {x$evalFail})))
return(list(pop=tmppop, fit=fit, evalFail=evalFail))
}

#' Configures the evaluation of the population of
#' a genetic algorithm.
#'
#' @description \code{xegaEvalPopulationFactory()} implements the selection
#'              of the evaluation function for the population of a 
#'              genetic algorithm.
#'
#' Current support:
#'
#' \enumerate{
#'   \item "EvalPopulation" returns 
#'        \code{xegaEvalPopulation}. (Default)
#'   \item "RepEvalPopulation" returns 
#'        \code{xegaReplEvalPopulation}. 
#'        For stochastic functions.
#'        Needs \code{lF$rep()} for the number of repetitions
#'        and \code{lF$apply()} for the (parallel) apply function.
#'    }
#'
#' @param method A string specifying the termination condition.
#'
#' @return A boolean function implementing the termination condition.
#'
#' @family Configuration
#'
#' @export
xegaEvalPopulationFactory<-function(method="EvalPopulation") {
if (method=="EvalPopulation") {f<-xegaEvalPopulation}
if (method=="RepEvalPopulation") {f<-xegaRepEvalPopulation}
if (!exists("f", inherits=FALSE))
        {stop("EvalPopulation label ", method, " does not exist")}
return(f)
}

