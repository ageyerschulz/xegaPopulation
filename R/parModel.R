#
# (c) 2021 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Population-level functions.
#                 Independent of gene representation.
#          Package: xegaPopulation.
#

# 
# Execution models of evaluation of a population:
# 

#' MultiCore apply of library parallel.
#'
#' @description 
#' The evaluation of the fitness of the genes of the population
#' is distributed to one worker on each core of the CPU of the 
#' local machine.
#' The package \code{parallel} of base R is used.
#' The number of cores is provided by \code{lF$Cores}. 
#' 
#' @details 
#' Be aware that
#' \itemize{
#' \item \code{parallel::mclapply()} assumes that each function evaluation 
#'       needs approximately the same time.
#' \item Best results are obtained if 
#'       \code{popsize} modulo \code{cores-1} is \code{0}.
#' \item Does not work on Windows.
#' }
#' 
#' @param pop        Population of genes.
#' @param EvalGene   Function for evaluating a gene.
#' @param lF         Local function configuration which provides 
#'                    all functions needed in \code{EvalGene()}.
#'
#' @return Fitness vector.
#' 
#' @family Execution Model
#' 
#' @examples
#' library(parallelly) 
#' if (supportsMulticore()){
#' lFxegaGaGene$Cores<-function() {2}
#' pop<-xegaInitPopulation(1000, lFxegaGaGene)
#' popnew<-MClapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @export
MClapply<-function(pop, EvalGene, lF) # nocov start
{       # z<-runif(1)
	parallel::mclapply(pop, 
		    EvalGene, 
		    lF=lF, 
		    mc.cores=max(1, lF$Cores()),
                    mc.set.seed = TRUE)
} # nocov end

#' MultiCore apply of library parallel for heterogenous tasks.
#'
#' @description 
#' The evaluation of the fitness of the genes of the population
#' is distributed to one worker on each core of the CPU of the 
#' local machine.
#' The package \code{parallel} of base R is used.
#' The number of cores is provided by \code{lF$Cores}. 
#' 
#' @details 
#' Be aware that
#' \itemize{
#' \item \code{parallel::mclapply()} assumes that each function evaluation 
#'       needs approximately the same time.
#' \item Best results are obtained if 
#'       \code{popsize} modulo \code{cores-1} is \code{0}.
#' \item Does not work on Windows.
#' }
#' 
#' @param pop        Population of genes.
#' @param EvalGene   Function for evaluating a gene.
#' @param lF         Local function configuration which provides 
#'                    all functions needed in \code{EvalGene()}.
#'
#' @return Fitness vector.
#' 
#' @family Execution Model
#' 
#' @examples
#' library(parallelly) 
#' if (supportsMulticore()){
#' lFxegaGaGene$Cores<-function() {2}
#' pop<-xegaInitPopulation(10, lFxegaGaGene)
#' popnew<-MClapplyHet(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @export
MClapplyHet<-function(pop, EvalGene, lF) # nocov start
{       # z<-runif(1)
	parallel::mclapply(pop, 
		    EvalGene, 
		    lF=lF, 
		    mc.cores=max(1, lF$Cores()),
                    mc.preschedule=FALSE,
                    mc.set.seed = TRUE)
} # nocov end

#' Future apply of R-package \code{future.apply}. 
#'
#' @description 
#' The \code{lapply()} function is redefined as as 
#' \code{future.apply::future_lapply()}.
#' Henrik Bengtsson recommends that the configuration of the 
#' parallel/distributed programming environment should be kept 
#' outside the package and left to the user. 
#' The advantage is that the user may take advantage of all 
#' parallel/distributed available backends for the Future API.
#' 
#' @details 
#' Be aware that
#' \itemize{
#' \item \code{future_lapply()} assumes that each function evaluation 
#'       need approximately the same time.
#' \item Best results are obtained 
#'       if \code{popsize} modulo \code{workers} is \code{0}.
#' }
#' 
#' @param pop        Population of genes.
#' @param EvalGene   Function for evaluating a gene.
#' @param lF         Local function factory which provides 
#'                    all functions needed in \code{EvalGene}.
#'
#' @return Fitness vector.
#' 
#' @references 
#' Bengtsson H (2021). “A Unifying Framework for Parallel and 
#' Distributed Processing in R using Futures.” 
#' The R Journal, 13(2), 208–227. <doi:10.32614/RJ-2021-048>
#'
#' @family Execution Model
#'
#' @examples 
#' pop<-xegaInitPopulation(1000, lFxegaGaGene)
#' library(future)
#' plan(multisession, workers=2)
#' popnew<-futureLapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' plan(sequential)
#'
#' @importFrom future.apply future_lapply
#' @export
futureLapply<-function(pop, EvalGene, lF) # nocov start
{
 future.apply::future_lapply(pop, EvalGene, lF=lF, future.seed=TRUE)
} # nocov end

#' Future apply of R-package \code{future.apply} configured
#' for a tasks with heterogenous execution times.
#'
#' @description 
#' The \code{lapply()} function is redefined as as 
#' \code{future.apply::future_lapply()}.
#' 
#' Henrik Bengtsson recommends that the configuration of the 
#' parallel/distributed programming environment should be kept 
#' outside the package and left to the user. 
#' The advantage is that the user may take advantage of all 
#' parallel/distributed available backends for the Future API.
#' 
#' @details 
#'   This configuration has an increased communication and 
#'   synchronization overhead.
#' 
#' @param pop        Population of genes.
#' @param EvalGene   Function for evaluating a gene.
#' @param lF         Local function factory which provides 
#'                    all functions needed in \code{EvalGene}.
#'
#' @return Fitness vector.
#' 
#' @references 
#' Bengtsson H (2021). “A Unifying Framework for Parallel and 
#' Distributed Processing in R using Futures.” 
#' The R Journal, 13(2), 208–227. <doi:10.32614/RJ-2021-048>
#'
#' @family Execution Model
#'
#' @examples 
#' pop<-xegaInitPopulation(30, lFxegaGaGene)
#' library(future)
#' plan(multisession, workers=2)
#' popnew<-futureLapplyHet(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' plan(sequential)
#'
#' @importFrom future.apply future_lapply
#' @export
futureLapplyHet<-function(pop, EvalGene, lF) # nocov start
{
 future.apply::future_lapply(pop, EvalGene, lF=lF, 
          future.seed=TRUE, future.chunk.size=NULL, future.scheduling=FALSE)
} # nocov end

#' uses parLapply of library parallel for using workers on 
#' machines in a local network. 
#'
#' @section Warning:
#'
#' This section has not been properly tested.
#' Random number generation?
#' Examples?
#'
#' @param pop        a population of genes.
#' @param EvalGene   the function for evaluating a gene.
#' @param lF          the local function factory which provides
#'                    all functions needed in \code{EvalGene}.
#'
#' @return Fitness vector.
#'
#' @family Execution Model
#'
#' @examples 
#' parm<-function(x) {function() {x}}
#' pop<-xegaInitPopulation(1000, lFxegaGaGene)
#' library(parallel)
#' clus<-makeCluster(2)
#' lFxegaGaGene$cluster<-parm(clus)
#' popnew<-PparLapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' stopCluster(clus)
#'
#' @importFrom parallel  clusterSetRNGStream
#' @importFrom parallel  parLapply
#' @export
PparLapply<-function(pop, EvalGene, lF) # nocov start
{       z<-runif(1)
        parallel::clusterSetRNGStream(lF$cluster(), NULL) # not reproducible.
	parallel::parLapply(lF$cluster(), pop, 
		    EvalGene, 
		    lF=lF)
} # nocov end

#### Low level specification of cluster object. Not recommended.
## clus<-makeCluster(spec=c("localhost", "localhost"), master="localhost", 
##                             port=1250, homogeneous=TRUE)

#' uses parLapplyLB of library parallel for using workers on 
#' machines in a local network. 
#'
#' @section Warning:
#'
#' This section has not been properly tested.
#' Random number generation?
#' Examples?
#'
#' @param pop        a population of genes.
#' @param EvalGene   the function for evaluating a gene.
#' @param lF          the local function factory which provides
#'                    all functions needed in \code{EvalGene}.
#'
#' @return Fitness vector.
#'
#' @family Execution Model
#'
#' @examples 
#' parm<-function(x) {function() {x}}
#' pop<-xegaInitPopulation(1000, lFxegaGaGene)
#' library(parallel)
#' clus<-makeCluster(2)
#' lFxegaGaGene$cluster<-parm(clus)
#' popnew<-PparLapplyHet(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#' stopCluster(clus)
#'
#' @importFrom parallel  clusterSetRNGStream
#' @importFrom parallel  parLapply
#' @export
PparLapplyHet<-function(pop, EvalGene, lF) # nocov start
{       z<-runif(1)
        parallel::clusterSetRNGStream(lF$cluster(), NULL) # not reproducible.
	parallel::parLapplyLB(lF$cluster(), pop, 
		    EvalGene, 
		    lF=lF)
} # nocov end
# to export data: see parallel::clusterExport

#' Configure the the execution model for gene evaluation.
#'
#' @description 
#' The current approach to distribution/parallelization of the genetic
#' algorithm is to parallelize the evaluation of the fitness function
#' only. The execution model defines the function \code{lF$lapply()}
#' used in the function \code{EvalPopulation()}.
#' 
#' @details
#' Currently we support the following parallelization models:
#' \enumerate{
#' \item "Sequential": Uses \code{base::lapply()}. (Default).
#' \item "MultiCore": Uses \code{parallel::mclapply()}. 
#'                    For tasks with approximately the same execution time.
#' \item "MultiCoreHet": Uses \code{parallel::mclapply()}. 
#'                    For tasks with a high variance of execution times.
#' \item "FutureApply": Uses \code{future.apply::future_lapply()}
#'                    Plans must be set up and
#'                    worker processes must be stopped.
#' \item "Cluster": Uses \code{parallel:parLapply()}.
#'                    A cluster object must be set up and the 
#'                    worker processes must be stopped. 
#' }
#'
#' The execution model \strong{"MultiCore"} provides parallelization restricted 
#' to a single computer: The master process starts R slave processes 
#' by fork() which are are run in separate memory spaces. 
#' At the time of fork() both memory spaces 
#' have the same content. Memory writes performed by one of the processes
#' do not affect the other. 
#'
#' The execution model \strong{"FutureApply"} makes the possibilities 
#' of the future backends for a wide range of parallel and distributed 
#' architectures available.  
#' The models of parallel resolving a future use 
#' different types of communication between master 
#' and slaves: 
#' \enumerate{ 
#' \item \code{plan(sequential)} configures sequential execution. Default.
#'
#' \item \code{w<-5; plan(multicore, workers=w)} configures an 
#'       asynchronous multicore execution of futures on 5 workers.
#'       
#' \item \code{w<-8; plan(multisession, workers=w)} configures a 
#'       multisession environment with 5 workers. 
#'       The evaluation of the future is done in parallel in 5 other 
#'       R sessions on the same machine. 
#'       Communication is done via socket connections, 
#'       the R sessions started serve multiple futures over their life time.
#'       The worker R sessions are stopped by calling \code{plan(sequential)}.
#'       The number of parallel sessions is restricted by the availability
#'       of connections. Up to R version 4.3, 
#'       a maximum of 125 connections is available.
#'
#' \item \code{w<-7; plan(callr, workers=w)} configures  
#'       the evaluation of futures on top of the \code{callr} package.
#'       The \code{callr} package creates for each future a separate R session.
#'       The communications is via files of serialized R objects.
#'       The advantages of \code{callr} are:
#'       \enumerate{  
#'       \item Each \code{callr} future is evaluated in a new R session
#'             which ends as soon as the value of the future has been 
#'             collected.
#'       \item The number of parallel \code{callr} futures is not restricted
#'             by the number of available connections, because the 
#'             communication is based on files of serialized R objects.
#'       \item No ports are used. This means no port clashes with other 
#'             processes and no firewall issues.
#'       }
#' 
#' \item Setting up a cluster environment for resolving futures works
#'       as follows. Write a function with the following elements:
#'       \enumerate{
#'       \item Generate a cluster object:
#'
#'        \code{cl<-makeClusterPSOCK(workers)} 
#'       \item Set up an on.exit condition for stopping the worker processes.
#'
#'       \code{on.exit(parallel::stopCluster(cl))}
#'       \item Set up the plan for resolving the future:
#'
#'       \code{oldplan<-plan(cluster, workers=cl)}
#'       \item Call the function with \code{future.apply::future_lapply}.
#'             E.g. the genetic algorithm.
#'       \item Restore the previous plan:
#'       \code{plan(oldplan)}
#'       }
#'       The cluster processes may be located on one or several computers.
#'       The communication between the processes is via sockets.
#'       Remote computers must allow the use of ssh to start R-processes
#'       without an interactive login.
#' }
#'
#' The execution model \strong{"Cluster"} allows the configuration of 
#' master-slave processing on local and remote machines.
#'
#' For evaluating tasks with highly variable execution times, 
#' it is recommended to use the corresponding heterogenous 
#' execution models which assign one task per computing node 
#' and start a new task to a node as soon as his task is finished.  
#' These execution models are "MultiCoreHet", "FutureApplyHet", and
#' "ClusterHet". Note that the communication and synchronization 
#' overhead of these execution models is substantially higher than 
#' for the homogenous execution models. 
#'
#' @param method   The label of the execution model: 
#'                  "Sequential" | 
#'                  "MultiCore" | "MultiCoreHet" |
#'                  "FutureApply" |  "FutureApplyHet" |
#'                  "Cluster" | "ClusterHet" .
#'
#' @return A function with the same result as the \code{lapply()}-function.  
#'
#' @family Configuration
#'
#' @export
ApplyFactory<- function(method="Sequential") {
if (method=="Sequential") {f<-base::lapply}
if (method=="MultiCore") {f<-MClapply}
if (method=="MultiCoreHet") {f<-MClapplyHet}
if (method=="FutureApply") {f<-futureLapply}
if (method=="FutureApplyHet") {f<-futureLapplyHet}
if (method=="Cluster") {f<-PparLapply}
if (method=="ClusterHet") {f<-PparLapplyHet}
if (!exists("f", inherits=FALSE))
        {stop("Execution Model label ", method, " does not exist")}
return(f)
}

