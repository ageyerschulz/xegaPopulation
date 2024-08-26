#
# (c) 2024 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-level functions.
#                 Independent of gene representation.
#          Package: xegaPopulation 
#

#' Exponential multiplicative cooling. 
#'
#' @description The temperature at time k is the net present value 
#'              of the starting temperature. The discount factor 
#'              is \code{lF$Alpha()}. 
#'              \code{lF$Alpha()} should be in \code{[0, 1]}.
#'
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$Alpha()} is the discount factor.
#'
#' @param k      Number of steps to discount.
#' @param lF     Local configuration.
#'
#' @return Temperature at time k.
#'
#' @references 
#'      Kirkpatrick, S., Gelatt, C. D. J, and Vecchi, M. P. (1983):
#'      Optimization by Simulated Annealing. 
#'      Science, 220(4598): 671-680.
#'      <doi:10.1126/science.220.4598.671>
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), Alpha=parm(0.99))
#' ExponentialMultiplicativeCooling(0, lF)
#' ExponentialMultiplicativeCooling(2, lF)
#' @export
ExponentialMultiplicativeCooling<-function(k, lF)
{ 
   return(lF$Temp0()*lF$Alpha()^k)	
}

#' Logarithmic multiplicative cooling. 
#'
#' @description This schedule decreases by the inverse proportion of the
#'              natural logarithm 
#'              of \code{k}. \code{lF$Alpha()} should be larger than 1.
#'
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$Alpha()} is a scaling factor.
#'
#' @param k      Number of steps to discount.
#' @param lF     Local configuration.
#'
#' @return Temperature at time k.
#'
#'      Aarts, E., and Korst, J. (1989):
#'      Simulated Annealing and Boltzmann Machines.
#'      A Stochastic Approach to Combinatorial Optimization and
#'      Neural Computing.
#'      John Wiley & Sons, Chichester.
#'      (ISBN:0-471-92146-7)
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), Alpha=parm(1.01))
#' LogarithmicMultiplicativeCooling(0, lF)
#' LogarithmicMultiplicativeCooling(2, lF)
#' @export
LogarithmicMultiplicativeCooling<-function(k, lF)
{ 
   return(lF$Temp0()/(1+lF$Alpha()*log(1+k)))	
}

#' Power multiplicative cooling. 
#'
#' @description This schedule decreases by the inverse proportion of 
#'              a power 
#'              of \code{k}. \code{lF$Alpha()} should be larger than 1.
#'
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          For \code{lF$CoolingPower()==1} and 
#'              \code{lF$CoolingPower()==2} this results in the 
#'          the linear and quadratic multiplicative cooling schemes 
#'          in A Comparison of Cooling Schedules for Simulated Annealing.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$Alpha()} is a scaling factor.
#'          \code{lF$CoolingPower()} is an exponential factor. 
#'
#' @param k    Number of steps to discount.
#' @param lF   Local configuration.
#'
#' @return Temperature at time k.
#'
#' @references The-Crankshaft Publishing (2023)
#'          A Comparison of Cooling Schedules for Simulated Annealing.
#'   <https://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/> 
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), Alpha=parm(1.01), CoolingPower=parm(2))
#' PowerMultiplicativeCooling(0, lF)
#' PowerMultiplicativeCooling(2, lF)
#' @export
PowerMultiplicativeCooling<-function(k, lF)
{ 
   return(lF$Temp0()/(1+lF$Alpha()*(k^lF$CoolingPower())))	
}

#' Power additive cooling. 
#'
#' @description This schedule decreases by a power of the
#'              n (= number of generations) linear fractions 
#'              between the starting temperature \code{lF$Temp0}
#'              and the final temperature \code{lF$tempN}. 
#'                   
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$TempN()} is the final temperature.
#'          \code{lF$CoolingPower()} is an exponential factor. 
#'          \code{lF$Generations()} is the number of generations (time). 
#'
#' @param k     Number of steps to discount.
#' @param lF    Local configuration.
#'
#' @return Temperature at time k.
#'
#' @references The-Crankshaft Publishing (2023)
#'          A Comparison of Cooling Schedules for Simulated Annealing.
#'   <https://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/> 
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), TempN=parm(10), Generations=parm(50), CoolingPower=parm(2))
#' PowerAdditiveCooling(0, lF)
#' PowerAdditiveCooling(2, lF)
#' @export
PowerAdditiveCooling<-function(k, lF)
{  fraction<-((lF$Generations()-k)/lF$Generations())^lF$CoolingPower()
   return(lF$TempN()+(lF$Temp0()-lF$TempN())*fraction)
}

#' Exponential additive cooling. 
#'
#' @description This schedule decreases in proportion to the
#'              inverse of \code{exp} raised to the 
#'              power of the temperature cycle in
#'              \code{lF$Generations()} (= number of generations) fractions 
#'              between the starting temperature \code{temp0}
#'              and the final temperature \code{tempN}. 
#'                   
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$TempN()} is the final temperature.
#'          \code{lF$Generations()} is the number of generations (time). 
#'
#' @param k      Number of steps to discount.
#' @param lF     Local configuration.
#'
#' @return The temperature at time k.
#'
#' @references The-Crankshaft Publishing (2023)
#'          A Comparison of Cooling Schedules for Simulated Annealing.
#'   <https://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/> 
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), TempN=parm(10), Generations=parm(50))
#' ExponentialAdditiveCooling(0, lF)
#' ExponentialAdditiveCooling(2, lF)
#' @export
ExponentialAdditiveCooling<-function(k, lF)
{  d<-lF$Temp0()-lF$TempN()
   t<-2*log(d)/lF$Generations()
   fraction<-1/(1+exp(t*(k-0.5*lF$Generations())))
   return(lF$TempN()+d*fraction)
}

#' Trigonometric additive cooling. 
#'
#' @description This schedule decreases in proportion to the
#'              cosine of the temperature cycle in
#'              \code{lF$Generations()} (= number of generations) fractions 
#'              between the starting temperature \code{lF$Temp0()}
#'              and the final temperature \code{lF$TempN()}. 
#'
#' @details Temperature is updated at the end of each generation
#'          in the main loop of the genetic algorithm.
#'          \code{lF$Temp0()} is the starting temperature.
#'          \code{lF$TempN()} is the final temperature.
#'          \code{lF$Generations()} is the number of generations (time). 
#'
#' @param k     Number of steps (time).
#' @param lF    Local configuration.
#'
#' @return Temperature at time k.
#'
#' @references The-Crankshaft Publishing (2023)
#'          A Comparison of Cooling Schedules for Simulated Annealing.
#'  <https://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/> 
#'
#' @family Cooling
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list(Temp0=parm(100), TempN=parm(10), Generations=parm(50))
#' TrigonometricAdditiveCooling(0, lF)
#' TrigonometricAdditiveCooling(2, lF)
#' @export
TrigonometricAdditiveCooling<-function(k, lF)
{  fraction<-0.5*(1+cos(k*pi/lF$Generations()))
   return(lF$TempN()+(lF$Temp0()-lF$TempN())*fraction)
}

#' Configure the cooling schedule of the acceptance 
#' function of a genetic algorithm.
#'
#' @description \code{CoolingFactory()} implements selection
#'              of a cooling schedule method.
#'
#' Current support:
#'
#' \enumerate{
#'   \item "ExponentialMultiplicative" returns 
#'        \code{ExponentialMultiplicativeCooling}. (Default)
#'   \item "LogarithmicMultiplicative" returns 
#'        \code{LogarithmicMultiplicativeCooling}.
#'   \item "PowerMultiplicative" returns 
#'        \code{PowerMultiplicativeCooling}.
#'        \code{coolingPower=1} specifies linear multiplicative cooling,  
#'        \code{coolingPower=2} specifies quadratic multiplicative cooling.
#'   \item "PowerAdditive" returns 
#'        \code{PowerAdditiveCooling}.
#'        \code{coolingPower=1} specifies linear additive cooling,  
#'        \code{coolingPower=2} specifies quadratic additive cooling.
#'   \item "ExponentialAdditive" returns 
#'        \code{ExponentialAdditiveCooling}.
#'   \item "TrigonometricAdditive" returns 
#'        \code{TrigonometricAdditiveCooling}.
#'    }
#'
#' @param method A string specifying the cooling schedule.
#'
#' @return A cooling schedule.
#'
#' @family Configuration
#'
#' @export
CoolingFactory<-function(method="ExponentialMultiplicative") {
if (method=="ExponentialMultiplicative") {f<- ExponentialMultiplicativeCooling}
if (method=="LogarithmicMultiplicative") {f<- LogarithmicMultiplicativeCooling}
if (method=="PowerMultiplicative") {f<- PowerMultiplicativeCooling}
if (method=="PowerAdditive") {f<- PowerAdditiveCooling}
if (method=="ExponentialAdditive") {f<- ExponentialAdditiveCooling}
if (method=="TrigonometricAdditive") {f<- TrigonometricAdditiveCooling}
if (!exists("f", inherits=FALSE))
        {stop("Cooling label ", method, " does not exist")}
return(f)
}

