library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

lFxegaGaGene$ReportEvalErrors<-function() {TRUE}

test_that("AcceptNewGene OK",
          {
	  id<-function(x, lF){x}
	  g1<-InitGene(lFxegaGaGene)
	  g2<-InitGene(lFxegaGaGene)
	  g2$evalFail<-TRUE
          expect_identical(AcceptNewGene(id, g1, lFxegaGaGene), g1)
          expect_identical(AcceptNewGene(id, g2, lFxegaGaGene), g2)
          }
)

test_that("AcceptBest OK",
          {
	  OPpipe1<-function(g, lF){InitGene(lF)}
	  set.seed(1)
	  g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
          lFfail<-lFxegaGaGene
          lFfail$penv$f<-function(parm, gene=0, lF=0) { "a"+sum(parm^{2})/0} 
          lFfail$CWorstFitness<-function() {-1000}
	  g2<-AcceptBest(OPpipe1, g1, lFfail)
	  g3<-AcceptBest(OPpipe1, g1, lFxegaGaGene)
	  g4<-AcceptBest(OPpipe1, g1, lFxegaGaGene)
	  g5<-AcceptBest(OPpipe1, g1, lFxegaGaGene)
	  g6<-AcceptBest(OPpipe1, g1, lFxegaGaGene)
          expect_identical(g1$fit, g2$fit)
          expect_identical(g1, g3)
          expect_identical(g1, g4)
          expect_identical(g1, g5)
          expect_identical(g6$fit>g1$fit, TRUE)
          }
)

test_that("AcceptMetropolis OK",
          {
          parm<-function(x){function() {return(x)}}
          lFxegaGaGene$Beta<-parm(1)
          lFxegaGaGene$TempK<-parm(10)
	  OPpipe1<-function(g, lF){InitGene(lF)}
	  set.seed(2)
	  g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
          lFfail<-lFxegaGaGene
          lFfail$penv$f<-function(parm, gene=0, lF=0) { "a"+sum(parm^{2})/0}
          lFfail$CWorstFitness<-function() {-1000}
          g2<-AcceptMetropolis(OPpipe1, g1, lFfail)
	  g3<-AcceptMetropolis(OPpipe1, g1, lFxegaGaGene)
	  g4<-AcceptMetropolis(OPpipe1, g1, lFxegaGaGene)
	  g5<-AcceptMetropolis(OPpipe1, g4, lFxegaGaGene)
          expect_identical(g1$fit, g2$fit)
          expect_identical(g1, g3)
          expect_identical((g1$fit>g4$fit), TRUE)
          expect_identical((g5$fit>g4$fit), TRUE)
          }
)

test_that("AcceptIVMetropolis OK",
          {
          parm<-function(x){function() {return(x)}}
          lFxegaGaGene$Beta<-parm(1)
          lFxegaGaGene$TempK<-parm(10)
          lFxegaGaGene$CBestFitness<-parm(20.39447)
          OPpipe1<-function(g, lF){InitGene(lF)}
          set.seed(2)
          g1<-lFxegaGaGene$EvalGene(InitGene(lFxegaGaGene), lFxegaGaGene)
          lFfail<-lFxegaGaGene
          lFfail$penv$f<-function(parm, gene=0, lF=0) { "a"+sum(parm^{2})/0}
          lFfail$CWorstFitness<-function() {-1000}
          g2<-AcceptIVMetropolis(OPpipe1, g1, lFfail)
          g3<-AcceptIVMetropolis(OPpipe1, g1, lFxegaGaGene)
          g4<-AcceptIVMetropolis(OPpipe1, g1, lFxegaGaGene)
          g5<-AcceptIVMetropolis(OPpipe1, g4, lFxegaGaGene)
          expect_identical(g1$fit, g2$fit)
          expect_identical(g1, g3)
          expect_identical((g1$fit>g4$fit), TRUE)
          expect_identical((g5$fit>g4$fit), TRUE)
          }
)

test_that("MetropolisAcceptanceProbability OK",
          {
          expect_equal(MetropolisAcceptanceProbability(0, 1, 20), 1)
          expect_equal(MetropolisAcceptanceProbability(1, 1, 20), 0.9512294,
	  tolerance=0.0001)
          expect_equal(MetropolisAcceptanceProbability(0, 5, 30), 1)
          expect_equal(MetropolisAcceptanceProbability(0, 0.5, 1000), 1)
          expect_equal(MetropolisAcceptanceProbability(1, 1, 20), 0.9512294,
	  tolerance=0.0001)
          expect_equal(MetropolisAcceptanceProbability(1, 2, 20), 0.9048374,
	  tolerance=0.0001)
          }
)

test_that("MetropolisTable OK",
          {
	  b<-c(0.6703200, 0.6676171, 0.6648980, 
	       0.6621626, 0.6594111, 0.6566433, 0.6538594,
               0.6510594, 0.6482432, 0.6454110)	       
	  a<-MetropolisTable(d=2, beta=2, T=10, alpha=0.99, steps=10)
          expect_equal(a[,5], b,
	  tolerance=0.0001)
          }
)

test_that("AcceptFactory OK",
          {
          Fun<-AcceptFactory()
          expect_identical(body(Fun), body(AcceptNewGene))
          Fun<-AcceptFactory("All")
          expect_identical(body(Fun), body(AcceptNewGene))
          Fun<-AcceptFactory("Best")
          expect_identical(body(Fun), body(AcceptBest))
          Fun<-AcceptFactory("Metropolis")
          expect_identical(body(Fun), body(AcceptMetropolis))
          Fun<-AcceptFactory("IVMetropolis")
          expect_identical(body(Fun), body(AcceptIVMetropolis))
          expect_error(AcceptFactory("Stchastic"))
          }
)

