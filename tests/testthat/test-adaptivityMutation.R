library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("ConstMRate OK",
	  {
	  parm<-function(x){function() {return(x)}}
	  lF<-list(MutationRate1=parm(0.20))
	  expect_identical(ConstMRate(100, lF), lF$MutationRate1())
	  }
)

test_that("IAMRate OK",
	  {
	  parm<-function(x){function() {return(x)}}
	  lF<-list(MutationRate1=parm(0.20), 
		   MutationRate2=parm(0.40), 
		   CutoffFit=parm(0.60), 
		   CBestFitness=parm(105))
	  expect_identical(IAMRate(100, lF), lF$MutationRate1())
	  expect_identical(IAMRate(50, lF), lF$MutationRate2())
	  }
)

test_that("MutationRateFactory OK",
          {
          Fun<-MutationRateFactory()
          expect_identical(body(Fun), body(ConstMRate))
          Fun<-MutationRateFactory("Const")
          expect_identical(body(Fun), body(ConstMRate))
          Fun<-MutationRateFactory("IV")
          expect_identical(body(Fun), body(IAMRate))
          expect_error(MutationRateFactory("Stchastic"))
          }
)

