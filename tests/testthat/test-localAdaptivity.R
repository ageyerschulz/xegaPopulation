library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("IAMBitRate OK",
	  {
	  parm<-function(x){function() {return(x)}}
	  lF<-list(BitMutationRate1=parm(0.20), 
		   BitMutationRate2=parm(0.40), 
		   CutoffFit=parm(0.60), 
		   CBestFitness=parm(105))
	  expect_identical(IAMBitRate(100, lF), lF$BitMutationRate1())
	  expect_identical(IAMBitRate(50, lF), lF$BitMutationRate2())
	  }
)

