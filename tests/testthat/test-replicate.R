library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("xegaNextPopulation P(Crossover=1.0), P(Mutation)=0.0 OK",
          {
          lFxegaGaGene$cGeneration<-parm(1)
          lFxegaGaGene$BitMutationRate1<-parm(0.001)
          lFxegaGaGene$MutationRate1<-parm(0.000)
	  lFxegaGaGene$MutationRate<-MutationRateFactory(method="Const")
	  lFxegaGaGene$ReplicateGene<-xegaGaReplicationFactory(method="Kid1")
	  lFxegaGaGene$CrossGene<-xegaGaCrossoverFactory(method="CrossGene")
	  lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
	  lFxegaGaGene$Accept<-AcceptFactory(method="All")
	   pop10<-xegaInitPopulation(10, lFxegaGaGene)
	   epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           npop10<-xegaNextPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
	   expect_equal(length(npop10), 10)
          }
)

test_that("xegaNextPopulation P(Crossover=0.0) P(Mutation=1.0)  OK",
          {
          lFxegaGaGene$cGeneration<-parm(1)
          lFxegaGaGene$BitMutationRate1<-parm(0.001)
          lFxegaGaGene$MutationRate1<-parm(1.0)
	  lFxegaGaGene$MutationRate<-MutationRateFactory(method="Const")
	  lFxegaGaGene$ReplicateGene<-xegaGaReplicationFactory(method="Kid1")
	  lFxegaGaGene$CrossGene<-xegaGaCrossoverFactory(method="CrossGene")
	  lFxegaGaGene$Elitist<-parm(FALSE)
	  lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
	  lFxegaGaGene$Accept<-AcceptFactory(method="All")
	   pop10<-xegaInitPopulation(10, lFxegaGaGene)
	   epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           npop10<-xegaNextPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
	   expect_equal(length(npop10), 10)
          }
)

test_that("xegaNextPopulation P(Crossover=0.0), P(Mutation)=0.0 OK",
          {
          lFxegaGaGene$cGeneration<-parm(1)
          lFxegaGaGene$BitMutationRate1<-parm(0.001)
          lFxegaGaGene$MutationRate1<-parm(0.0)
	  lFxegaGaGene$MutationRate<-MutationRateFactory(method="Const")
	  lFxegaGaGene$ReplicateGene<-xegaGaReplicationFactory(method="Kid1")
	  lFxegaGaGene$CrossGene<-xegaGaCrossoverFactory(method="CrossGene")
	  lFxegaGaGene$CrossRate<-function(fit, lF) {0.0}
	  lFxegaGaGene$Accept<-AcceptFactory(method="All")
	   pop10<-xegaInitPopulation(10, lFxegaGaGene)
	   epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           npop10<-xegaNextPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
	   expect_equal(length(npop10), 10)
          }
)

test_that("xegaNextPopulation P(Crossover=1.0), P(Mutation)=1.0 OK",
          {
          lFxegaGaGene$cGeneration<-parm(1)
          lFxegaGaGene$BitMutationRate1<-parm(0.001)
          lFxegaGaGene$MutationRate1<-parm(1.0)
	  lFxegaGaGene$MutationRate<-MutationRateFactory(method="Const")
	  lFxegaGaGene$ReplicateGene<-xegaGaReplicationFactory(method="Kid1")
	  lFxegaGaGene$CrossGene<-xegaGaCrossoverFactory(method="CrossGene")
	  lFxegaGaGene$CrossRate<-function(fit, lF) {1.0}
	  lFxegaGaGene$Accept<-AcceptFactory(method="All")
	   pop10<-xegaInitPopulation(10, lFxegaGaGene)
	   epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           npop10<-xegaNextPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
	   expect_equal(length(npop10), 10)
          }
)

