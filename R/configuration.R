#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Population-level functions.
#                 Independent of gene representation.
#          Package: xegaPopulation.
#

#' Remembers R command command with which algorithm has been called.
#' 
#' @description \code{xegaConfiguration()} returns the command with which
#'              the genetic algorithm has been called. 
#'              For replicating computational experiments with 
#'              genetic algorithms.
#' 
#' @param GAname    Name of genetic algorithm's main function.
#'                      (Currently: "Run").
#' @param penv      The expression for the problem environment \code{penv}.
#'                       Use: \code{substitute(penv)}.
#' @param grammar   The grammar \code{grammar}.
#'                       Use: \code{substitute(grammar)}.
#' @param env       Environment with variable value bindings.
#'                       Use: \code{environment()}.
#'
#' @return A named list with the following elements:
#'   \itemize{
#'   \item \code{$GAconf}: 
#'                  A text string with the call of the genetic algorithm
#'                  (the function we want to capture the call).
#'   \item \code{$GAenv}:   The environment with the arguments bound to the 
#'                  values when the genetic algorithm was called.
#'   }
#' 
#' @section Warning:
#'    \itemize{
#'    \item
#'    $GAenv is correct only for simple arguments (strings or numbers) 
#'    not for complex objects like problem environments.
#'    \item 
#'    \code{future.apply::future_lapply()} is configured by a plan 
#'    statement which must be issued before calling the genetic 
#'    algorithm. At the moment, the plan chosen is not remembered.
#'    } 
#'
#' @family Configuration
#'
#' @examples
#' GA<-function(pe, grammar=NULL, nope=1.5, sle="test", ok=TRUE) 
#' {xegaConfiguration("GA", substitute(pe), substitute(grammar), environment())}
#' Para<-5
#' GA(Para)
#' Cube<-7
#' GA(Cube, 2, 3, 4)
#'
#' @export
xegaConfiguration<-function(GAname, penv, grammar, env)
{       
arglst<-as.list(env)
penvName<-deparse(penv)
if (is.null(grammar))
{grammarName<-"NULL"} else
{grammarName<-deparse(grammar)}
GAcall<-paste(GAname,"(penv=", penvName, ",grammar=", grammarName, sep="")
for (i in 3: length(arglst))
{
if (typeof(unlist(arglst[i]))=="character")
{GAcall<-paste(GAcall, ",", names(arglst[i]), "=\"", arglst[i], "\"", sep="")}
else
{GAcall<-paste(GAcall, ",", names(arglst[i]), "=", arglst[i], sep="")}
}       
GAcall<-paste(GAcall,")", sep="")

return(list(GAconf=GAcall, GAenv=arglst))
}  
