library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("ConstCRate OK",
	  {
	  parm<-function(x){function() {return(x)}}
	  lF<-list(CrossRate1=parm(0.20))
	  expect_identical(ConstCRate(100, lF), lF$CrossRate1())
	  }
)

test_that("IACRate OK",
	  {
	  parm<-function(x){function() {return(x)}}
	  lF<-list(CrossRate1=parm(0.20), 
		   CrossRate2=parm(0.40), 
		   CutoffFit=parm(0.60), 
		   CBestFitness=parm(105))
	  expect_identical(IACRate(100, lF), lF$CrossRate1())
	  expect_identical(IACRate(50, lF), lF$CrossRate2())
	  }
)

test_that("CrossRateFactory OK",
          {
          Fun<-CrossRateFactory()
          expect_identical(body(Fun), body(ConstCRate))
          Fun<-CrossRateFactory("Const")
          expect_identical(body(Fun), body(ConstCRate))
          Fun<-CrossRateFactory("IV")
          expect_identical(body(Fun), body(IACRate))
          expect_error(CrossRateFactory("Stchastic"))
          }
)

