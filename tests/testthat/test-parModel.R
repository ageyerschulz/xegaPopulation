library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

#test_that("MClapply OK",
#  {
#        skip_on_cran()
#	pop<-xegaInitPopulation(1000, lFxegaGaGene)
#        npop0<-lapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#        npop1<-MClapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
#        expect_identical(npop0, npop1)
#  })

test_that("ApplyFactory OK",
          {
          Fun<-ApplyFactory()
          expect_identical(body(Fun), body(base::lapply))
          Fun<-ApplyFactory("Sequential")
          expect_identical(body(Fun), body(base::lapply))
          Fun<-ApplyFactory("MultiCore")
          expect_identical(body(Fun), body(MClapply))
          #Fun<-ApplyFactory("FutureMultiCore")
          #expect_identical(body(Fun), body(lapplyMCFuture))
          Fun<-ApplyFactory("Cluster")
          expect_identical(body(Fun), body(PparLapply))
          expect_error(ApplyFactory("Stchastic"))
          }
)

