library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("terminatedFalse OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
                   Alpha=parm(0.99))
          expect_identical(terminatedFalse(0, lF), FALSE)
          }
)

test_that("TerminationFactory OK",
          {
          Fun<-TerminationFactory()
          expect_identical(body(Fun), body(terminatedFalse))
          Fun<-TerminationFactory("NoTermination")
          expect_identical(body(Fun), body(terminatedFalse))
          expect_error(TerminationFactory("Stchastic"))
          }
)

