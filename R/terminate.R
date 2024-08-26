

#
# (c) 2024 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package: xegaPopulation
#

#' No termination condition.
#' 
#' @description A boolean function which always returns \code{FALSE}.
#'
#' @param solution  A named list with at least the following elements:
#'                  $name, $fitness, $value, $numberOfSolutions, 
#'                  $genotype, $phenotype, $phenotypeValue.
#'                  
#' @param lF        Local function configuration.                   
#'
#' @return  FALSE
#' 
#' @family Termination Condition
#'
#' @examples
#'    lF<-list()
#'    terminatedFalse(1.0, lF)
#' @export
terminatedFalse<-function(solution, lF) {FALSE}

#' Check terminatedFalse()
#' 
#' @param penv    A problem environment.
#' @param max     Maximize?
#'
#' @return A named list
#'         \itemize{
#'          \item \code{$OK}   \code{TRUE}
#'          \item \code{$penv} \code{penv}
#'                 }
#'
checkTerminatedFalse<-function(penv, max) 
   {lst<-list(); lst$OK<-TRUE; lst$penv<-penv; return(lst)}

# For benchmarks with known global optima:

#' Terminates, if the absolute deviation from the global optimum is small.
#'
#' @description \code{terminateAbsoluteError()} 
#' returns \code{TRUE} if the value of the current solution 
#' is in the interval from (globalOptimum - eps) to 
#' (globalOptimum + eps). 
#' 
#' @details Useful for benchmark functions with known global optima.
#'          
#'
#' @param solution  A named list with at least the following elements:
#'                  $name, $fitness, $value, $numberOfSolutions, 
#'                  $genotype, $phenotype, $phenotypeValue.
#'                  
#' @param lF        Local function configuration. It must contain
#'                  \itemize{
#'                  \item \code{lF$penv$globalOptimum()} which returns 
#'                        the global optimum. 
#'                  \item \code{lF$TerminationEps()} which specifies the 
#'                        the maximal allowed deviation of the current 
#'                        best solution from the global optimum.
#'                  }                   
#'
#' @return Boolean.
#' 
#' @family Termination Condition
#'
#' @examples
#'     parm<-function(x){function() {return(x)}}
#'     olst<-list(); olst$value<-10
#'     penv<-list(); penv$globalOptimum<-parm(olst)
#'     lF<-list(); lF$penv<-penv; lF$TerminationEps<-parm(1.2);lF$Max<-parm(1.0)
#'     solution<-list(); solution$genotype<-list(); solution$genotype$fit<-8.0
#'     terminateAbsoluteError(solution, lF)
#'     solution<-list(); solution$genotype<-list(); solution$genotype$fit<-8.9
#'     terminateAbsoluteError(solution, lF)
#' @export
terminateAbsoluteError<-function(solution, lF) {
             opt<-lF$penv$globalOptimum()$value
             eps<-abs(lF$TerminationEps())
            # if ((solution$phenotypeValue>(opt-eps)) &
            #    (solution$phenotypeValue<(opt+eps)))
             copt<-solution$genotype$fit*lF$Max()
#      cat("copt:", copt, "opt-eps", (opt-eps), "opt+eps", (opt+eps), "\n")
             if ((copt>(opt-eps)) &
                (copt<(opt+eps)))
             {return(TRUE)} else {return(FALSE)}
}

#' Check terminateError()
#' 
#' @param penv    A problem environment.
#' @param max     Maximize?
#'
#' @return A named list
#'         \itemize{
#'          \item \code{$OK}   \code{TRUE}
#'          \item \code{$penv} \code{penv}
#'                 }
#'
checkTerminateError<-function(penv, max) 
{ nms<-names(penv)
  parm<-function(x){function() {return(x)}}
     if ("globalOptimum" %in% nms) 
        {lst<-list(); lst$OK<-TRUE; lst$penv<-penv}
     if ((!("globalOptimum" %in% nms)) 
        && ("solution" %in% nms))
        {
         npenv<-penv 
         # print(names(penv$solution()))
         # print(("maximum" %in% names(penv$solution())))
         if (max==TRUE) 
         {
            if ("maximum" %in% names(penv$solution()))
            { olst<-list(); olst$value<-penv$solution()$maximum
              npenv$globalOptimum<-parm(olst)
              lst<-list(); lst$OK<-TRUE; lst$penv<-npenv}
            else
            {stop("maximum not specified in penv$solution().")}
         }
         if (max==FALSE)
         {
            if ("minimum" %in% names(penv$solution()))
            { olst<-list(); olst$value<-penv$solution()$minimum
              npenv$globalOptimum<-parm(olst)
              lst<-list(); lst$OK<-TRUE; lst$penv<-npenv}
            else
            {stop("minimum not specified in penv$solution().")}
         }
        }
     if ((!("globalOptimum" %in% nms)) 
        && (!("solution" %in% nms)))
        { a<-"Global optimum not known.\n"
          b<-"penv must provide the global optimum value\n" 
          c<-"either in penv$globalOptimum() or in penv$solution()!"
          stop(a, b, c)}
     return(lst)
}          

#' Terminates, if the relative deviation from the global optimum is small.
#'
#' @description \code{terminateRelativeError()} 
#' returns \code{TRUE} if the value of the current solution 
#' is in the interval from (globalOptimum - (globalOptimum*eps)) to 
#' (globalOptimum + (globalOptimum*eps)). 
#' 
#' @details Useful for benchmark functions with known global optima.
#'          Note that for a global optimum of \code{0} this function 
#'          fails.
#'          
#'
#' @param solution  A named list with at least the following elements:
#'                  $name, $fitness, $value, $numberOfSolutions, 
#'                  $genotype, $phenotype, $phenotypeValue.
#'                  
#' @param lF        Local function configuration. It must contain
#'                  \itemize{
#'                  \item \code{lF$penv$globalOptimum()} which returns 
#'                        the global optimum. 
#'                  \item \code{lF$TerminationEps()} which specifies the 
#'                        the fraction of the global optimum
#'                        used for computing the upper and lower bounds
#'                        for the interval in which the best current 
#'                        solution must be for terminating the algorithm. 
#'                  }                   
#'
#' @return Boolean.
#' 
#' @family Termination Condition
#'
#' @examples
#'     parm<-function(x){function() {return(x)}}
#'     olst<-list(); olst$value<-10
#'     penv<-list(); penv$globalOptimum<-parm(olst)
#'     lF<-list(); lF$penv<-penv; lF$TerminationEps<-parm(1.2);lF$Max<-parm(1.0)
#'     solution<-list(); solution$genotype<-list();  solution$genotype$fit<-8.0
#'     terminateRelativeError(solution, lF)
#'     solution<-list(); solution$genotype<-list();  solution$genotype$fit<-9.6
#'     terminateRelativeError(solution, lF)
#' @export
terminateRelativeError<-function(solution, lF) {
             opt<-lF$penv$globalOptimum()$value
             eps<-abs(lF$TerminationEps()*opt)
             copt<-solution$genotype$fit*lF$Max()
             if ((copt>(opt-eps)) &
                (copt<(opt+eps)))
             # if ((solution$phenotypeValue>(opt-eps)) &
             #   (solution$phenotypeValue<(opt+eps)))
             {return(TRUE)} else {return(FALSE)}
}

#' Terminates if relative deviation from optimum is small. Works at 0.
#'
#' @description \code{terminateRelativeErrorZero()} 
#' returns \code{TRUE} if the value of the current solution 
#' is in the interval from (globalOptimum - (globalOptimum*eps)) to 
#' (globalOptimum + (globalOptimum*eps)). If globalOptimum is zero, 
#' test interval (0-eps) to (0+eps).
#' 
#' @details Useful for benchmark functions with known global optima.
#'          Note that for a global optimum of \code{0} this function 
#'          terminates if the current optimum is between 
#'          \code{0-terminationEps} and \code{0+terminationEps}. 
#'          
#' @param solution  A named list with at least the following elements:
#'                  $name, $fitness, $value, $numberOfSolutions, 
#'                  $genotype, $phenotype, $phenotypeValue.
#'                  
#' @param lF        Local function configuration. It must contain
#'                  \itemize{
#'                  \item \code{lF$penv$globalOptimum()} which returns 
#'                        the global optimum. 
#'                  \item \code{lF$TerminationEps()} which specifies the 
#'                        the fraction of the global optimum
#'                        used for computing the upper and lower bounds
#'                        for the interval in which the best current 
#'                        solution must be for terminating the algorithm. 
#'                  }                   
#'
#' @return Boolean.
#' 
#' @family Termination Condition
#'
#' @examples
#'     parm<-function(x){function() {return(x)}}
#'     olst<-list(); olst$value<-0
#'     penv<-list(); penv$globalOptimum<-parm(olst)
#'     lF<-list(); lF$penv<-penv; lF$TerminationEps<-parm(1.2);lF$Max<-parm(1.0)
#'     solution<-list(); solution$genotype<-list();  solution$genotype$fit<-0.5
#'     terminateRelativeErrorZero(solution, lF)
#'     solution<-list(); solution$genotype<-list();  solution$genotype$fit<-9.6
#'     terminateRelativeErrorZero(solution, lF)
#' @export
terminateRelativeErrorZero<-function(solution, lF) {
             opt<-lF$penv$globalOptimum()$value
             eps<-max(abs(lF$TerminationEps()*opt), abs(lF$TerminationEps()))
             copt<-solution$genotype$fit*lF$Max()
             if ((copt>(opt-eps)) &
                (copt<(opt+eps)))
             # if ((solution$phenotypeValue>(opt-eps)) &
             #   (solution$phenotypeValue<(opt+eps)))
             {return(TRUE)} else {return(FALSE)}
}

#' Configure the termination condition(s) 
#' a genetic algorithm.
#'
#' @description \code{TerminationFactory()} implements the selection
#'              of a termination method.
#'
#' Current support:
#'
#' \enumerate{
#'   \item "NoTermination" returns 
#'        \code{terminatedFalse}. (Default)
#'   \item "AbsoluteError" returns 
#'        \code{terminateAbsoluteError()}. 
#'        For benchmark functions with known global optima.
#'        Termination condition is fulfilled if the current 
#'        best solution is in the interval 
#'        from (globalOptimum-eps) to (globalOptimum+eps).
#'   \item "RelativeError" returns 
#'        \code{terminateRelativeError()}. 
#'        For benchmark functions with known global optima.
#'        Termination condition is fulfilled if the current 
#'        best solution is in the interval 
#'        from (globalOptimum-(globalOptimum*eps)) to 
#'        (globalOptimum+(globalOptimum*eps)).
#'        Does not specify an interval if globalOptimum is zero.
#'   \item "RelativeErrorZero" returns 
#'        \code{terminateRelativeErrorZero()}. 
#'        For benchmark functions with known global optima.
#'        Termination condition is fulfilled if the current 
#'        best solution is in the interval 
#'        from (globalOptimum-(globalOptimum*eps)) to 
#'        (globalOptimum+(globalOptimum*eps)).
#'        If the globalOptimum is zero, the interval is 
#'        from \code{-terminationEps} to \code{terminationEps}. 
#'    }
#'
#' @param method A string specifying the termination condition.
#'
#' @return A boolean function implementing the termination condition.
#'
#' @family Configuration
#'
#' @export
TerminationFactory<-function(method="NoTermination") {
if (method=="NoTermination") {f<-terminatedFalse}
if (method=="AbsoluteError") {f<-terminateAbsoluteError}
if (method=="RelativeError") {f<-terminateRelativeError}
if (method=="RelativeErrorZero") {f<-terminateRelativeErrorZero}
if (!exists("f", inherits=FALSE))
        {stop("Termination label ", method, " does not exist")}
return(f)
}

#' Configure consistency checks and adapt \code{penv} for terminationConditions.
#'
#' @description For each termination condition, a check must be provided.
#'              A check fails (stops) if the consistency requirements 
#'              of a termination condition are not fulfilled.
#'              However, a check may modify the problem environment to 
#'              establish consistency.
#'
#' @param method A string specifying the termination condition.
#'
#' @return A check function. 
#'
#' @family Configuration
#'
#' @export
checkTerminationFactory<-function(method="NoTermination") {
if (method=="NoTermination") {f<-checkTerminatedFalse}
if (method=="AbsoluteError") {f<-checkTerminateError}
if (method=="RelativeError") {f<-checkTerminateError}
if (method=="RelativeErrorZero") {f<-checkTerminateError}
if (!exists("f", inherits=FALSE))
        {stop("Termination label ", method, " does not exist")}
return(f)
}
