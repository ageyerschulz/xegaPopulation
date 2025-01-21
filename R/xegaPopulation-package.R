
#' Population level functions 
#' 
#' The \code{xegaPopulation} package provides the representation independent
#' functions of the population level of the simple genetic 
#' algorithm xegaX packages: 
#' \itemize{
#' \item File xegaPopulation.R:
#' \itemize{
#' \item Initializing a population of genes.
#' \item Getting the indices of the best genes in a population of genes
#'       for getting the best solution(s) in a population of genes.
#' \item Configurable summary report of population fitness statistics.
#' \item Observation of the summary statistics of a population of genes.
#' \item Logging of the phenotype and the value of the phenotype.
#' }
#' \item File xegaNextPopulation.R:
#' \itemize{
#' \item Computation of the next population of genes.
#' \item Evaluation of the next population of genes.
#' }
#'
#' \strong{Future}: Improved support for parallelization suggests a 
#'                  different division of labor: 
#'   \itemize{
#'       \item Construct a list of abstract task descriptions 
#'             with one element per gene.
#'       \item Provide for a parallel execution of these task descriptions. 
#'             This requires changes in the structuring of the 
#'             operator pipelines and the replicate gene functions 
#'             for the different gene representations and algorithms.
#'       \item Performance improvement depends on the gene representation 
#'             and on the use of function evaluations in the genetic 
#'             machinery. For example, for the TSP problem, 
#'             function evaluations 
#'             are embedded into most of the mutation operators.
#'           }
#'
#' \item File acceptance.R: 
#'       Acceptance rules for new genes and a function factory for configuring
#'       them. 
#' \item File cooling.R: Cooling schedules for temperature reduction.
#'
#' \item File localAdaptivity.R: Unused. 
#'        Move to gene dependent packages planned.
#' \item File adaptivityCrossover.R: 
#'               Functions constant and adaptive crossover rates. 
#' \item File adaptivityMutation.R:
#'               Functions constant and adaptive mutation rates. 
#' \item File parModel.R: Execution models for parallelization.
#' \itemize{
#' \item "Sequential": Configures lapply as \code{lapply()}.
#' \item "MultiCore": Configures lapply as \code{parallel::mclapply()}.
#'                    The number of cores is set by \code{lF$Core()}. 
#' }
#' \item File configuration.R: Documenting how the algorithm was called.
#'                         Support for the replication of computational
#'                         experiments (replicate and replay).
#' }
#'
#' @section Interface of Acceptance Rules:
#' 
#' \code{newGene<-accept(OperatorPipeline, gene, lF)}
#'
#' \enumerate{
#' \item Accept all new genes: Identity function. For genetic algorithms.
#' \item Accept best: Accepts the gene with the highest fitness.
#'                       For greedy and randomized greedy algorithms
#'                       (hill-climbing algorithms).
#' \item The Metropolis and the individually variable Metropolis rule:
#'       If the new gene gene is better, accept it.
#'       If the old gene is better, make a biased random choice. 
#'       The probability of accepting a decrease in fitness depends on 
#'       the fitness distance between genes, a constant beta for scaling 
#'       the exponential decay and a temperature parameter and for 
#'       the individually variable Metropolis rule a correction term 
#'       which depends on the distance to the best known fitness of the run. 
#' }
#'
#' \strong{Constants for Acceptance Rules.}
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$Beta()         \tab  ?    \tab AcceptMetropolis() \cr 
#'                   \tab       \tab AcceptIVMetropolis() \cr 
#' lF$TempK()        \tab  ?    \tab AcceptMetropolis() \cr 
#'                   \tab       \tab AcceptIVMetropolis() \cr 
#' lF$lFCBestFitness() \tab None \tab AcceptIVMetropolis() \cr 
#' }
#'
#' @section Interface of Cooling Schedules:
#'
#' \code{Temperature<-cooling(k, lF)}
#'
#' Cooling schedules convert the progress of the time in the algorithm
#' (measured in generations) into a temperature.
#' The temperature influences the probability of accepting a gene
#' with less fitness than its parent gene.
#'
#' \strong{Constants for Cooling Schedules.}
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$Alpha()         \tab  ?    \tab ExponentialMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab LogarithmicMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab PowerMultiplicativeCooling() \cr 
#' lF$Temp0()         \tab  ?    \tab ExponentialMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab LogarithmicMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab PowerMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab PowerAdditiveCooling() \cr 
#'                    \tab  ?    \tab ExponentialAdditiveCooling() \cr 
#'                    \tab  ?    \tab TrigonometricAdditiveCooling() \cr 
#' lF$TempN()         \tab  ?    \tab PowerAdditiveCooling() \cr 
#'                    \tab  ?    \tab ExponentialAdditiveCooling() \cr 
#'                    \tab  ?    \tab TrigonometricAdditiveCooling() \cr 
#' lF$CoolingPower()  \tab  ?    \tab PowerMultiplicativeCooling() \cr 
#'                    \tab  ?    \tab PowerAdditiveCooling() \cr 
#' lF$Generations()   \tab       \tab PowerAdditiveCooling() \cr 
#'                    \tab       \tab ExponentialAdditiveCooling() \cr 
#'                    \tab  ?    \tab TrigonometricAdditiveCooling() \cr 
#' }
#'
#' @section Interface of Rates:
#'
#' \code{rate<-rateFunction(fit, lF)}
#'
#' Crossover and mutation rate functions may be adaptive.
#' The interface allows for dependencies of the rate 
#' on fitness and constants in the local configuration.
#'
#' \strong{Constants for Adaptive Crossover and Mutation Rates}
#'
#' \tabular{rcl}{ 
#' \strong{Constant}  \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$CrossRate1()    \tab  ?    \tab IACRate() \cr 
#' lF$CrossRate2()    \tab  ?    \tab IACRate() \cr 
#' lF$MutationRate1() \tab       \tab IAMRate() \cr 
#' lF$MutationRate2() \tab       \tab IAMRate() \cr 
#' lF$CutoffFit()     \tab  ?    \tab IACRate() \cr 
#' lF$CBestFitness()  \tab       \tab IACRate() \cr 
#'                    \tab       \tab IAMRate() \cr 
#' }
#'
#' @section Interface of Termination Conditions:
#'
#' \code{hasTerminated<-Terminate(solution, lF)}
#'
#' The interface allows the specification of termination conditions 
#' for the genetic algorithm. The abstract function \code{Terminate}
#' returns a boolean value. TBD
#' 
#' \strong{Dependencies of Termination Conditions}
#' 
#' \tabular{rr}{ 
#' \strong{Condition}  \tab \strong{Requires} \cr 
#' terminatedFalse()   \tab -                 \cr
#' terminateAbsoluteError() \tab lF$penv$globalOptimum() \cr
#'                          \tab lF$TerminationEps \cr
#' terminateRelativeError() \tab lF$penv$globalOptimum() \cr
#'                          \tab lF$TerminationEps \cr
#' }
#'
#' @section The Architecture of the xegaX-Packages:
#' 
#' The xegaX-packages are a family of R-packages which implement 
#' eXtended Evolutionary and Genetic Algorithms (xega).  
#' The architecture has 3 layers, 
#' namely the user interface layer,
#' the population layer, and the gene layer: 
#' 
#' \itemize{
#' \item
#' The user interface layer (package \code{xega}) 
#' provides a function call interface and configuration support
#' for several algorithms: genetic algorithms (sga), 
#' permutation-based genetic algorithms (sgPerm), 
#' derivation free algorithms as e.g. differential evolution (sgde), 
#' grammar-based genetic programming (sgp) and grammatical evolution
#' (sge). 
#'
#' \item
#' The population layer (package \code{xegaPopulation}) contains
#' population related functionality as well as support for 
#' population statistics dependent adaptive mechanisms and parallelization.
#'
#' \item 
#' The gene layer is split in a representation independent and 
#' a representation dependent part:
#' \enumerate{
#' \item 
#'  The representation indendent part (package \code{xegaSelectGene})
#'  is responsible for variants of selection operators, evaluation 
#'  strategies for genes, as well as profiling and timing capabilities.        
#' \item 
#'  The representation dependent part consists of the following packages: 
#' \itemize{
#' \item \code{xegaGaGene} for binary coded genetic algorithms.
#' \item \code{xegaPermGene} for permutation-based genetic algorithms.
#' \item \code{xegaDfGene} for derivation free algorithms as e.g. 
#'                         differential evolution.
#' \item \code{xegaGpGene} for grammar-based genetic algorithms.
#' \item \code{xegaGeGene} for grammatical evolution algorithms.
#' }
#' The packages \code{xegaDerivationTrees} and \code{xegaBNF} support
#' the last two packages:
#' \code{xegaBNF} essentially provides a grammar compiler and 
#' \code{xegaDerivationTrees} an abstract data type for derivation trees.
#' }} 
#' 
#' @family Package Description 
#'
#' @name xegaPopulation
#' @aliases xegaPopulation
#' @docType package
#' @title Package xegaPopulation.
#' @author Andreas Geyer-Schulz
#' @section Copyright: (c) 2023 Andreas Geyer-Schulz
#' @section License: MIT
#' @section URL: https://github.com/ageyerschulz/xegaPopulation
#' @section Installation: From CRAN by \code{install.packages('xegaPopulation')}
"_PACKAGE"
