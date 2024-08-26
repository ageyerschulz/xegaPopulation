library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("lFxegaGaGene OK",
          {
           lFxegaGaGene$CrossRate<-parm(0.2)
	   expect_identical(lFxegaGaGene$penv$name(), "Parabola2D")
           expect_equal(lFxegaGaGene$replay(), 0)
           expect_equal(lFxegaGaGene$verbose(), 4)
           expect_equal(lFxegaGaGene$CutoffFit(), 0.5)
           expect_equal(lFxegaGaGene$CBestFitness(), 100)
           expect_equal(lFxegaGaGene$MutationRate1(), 0.01)
           expect_equal(lFxegaGaGene$MutationRate2(), 0.20)
           expect_equal(lFxegaGaGene$CrossRate(), 0.2)
           expect_equal(lFxegaGaGene$UCrossSwap(), 0.2)
           expect_equal(lFxegaGaGene$Max(), 1)
           expect_equal(lFxegaGaGene$Offset(), 1)
           expect_equal(lFxegaGaGene$Eps(), 0.01)
           expect_identical(lFxegaGaGene$Elitist(), TRUE)
           expect_equal(lFxegaGaGene$TournamentSize(), 2)
          }
)

test_that("xegaInitPopulation OK",
          {
	   pop10<-xegaInitPopulation(10, lFxegaGaGene)
           expect_equal(length(pop10), 10)
 	   fit10<-unlist(lapply(pop10, function(x) {x$fit}))
           expect_identical(fit10, rep(0,10))
 	   ev10<-unlist(lapply(pop10, function(x) {x$evaluated}))
           expect_identical(ev10, rep(FALSE,10))
          }
)

test_that("xegaObservePopulation (Start population OK",
          {
           pop10<-xegaInitPopulation(10, lFxegaGaGene)
           fit10<-unlist(lapply(pop10, function(x) {x$fit}))
           expect_identical(xegaObservePopulation(fit10), rep(0,8))
          }
)

test_that("xegaObservePopulation (Evaluated population OK",
          {
           pop10<-xegaInitPopulation(10, lFxegaGaGene)
           epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           expect_equal(xegaObservePopulation(epop10$fit)[1], mean(epop10$fit))
          }
)

test_that("xegaEvalPopulation OK",
          {
           pop10<-xegaInitPopulation(10, lFxegaGaGene)
           epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           expect_identical(names(epop10), c("pop", "fit", "evalFail"))
           expect_equal(length(epop10$pop), 10)
           expect_equal(length(epop10$fit), 10)
          }
)

### Not on cran!
test_that("xegaEvalPopulation MultiCore OK",
          {
	   skip_on_cran()
	   lFxegaGaGene[["evalPopLapply"]]<-ApplyFactory(method="MultiCore")
           pop10<-xegaInitPopulation(10, lFxegaGaGene)
           epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
           expect_identical(names(epop10), c("pop", "fit", "evalFail"))
           expect_equal(length(epop10$pop), 10)
           expect_equal(length(epop10$fit), 10)
          }
)

test_that("xegaBestGeneInPopulation 1 max OK",
          {
           f<-sample(1000, 100)
	   index<-which(max(f)==f)
           expect_identical(xegaBestGeneInPopulation(f), index)
          }
)

test_that("xegaBestGeneInPopulation 1 max OK",
          {
           f<-sample(20, 100, replace=TRUE)
	   index<-which(max(f)==f)
           expect_identical(xegaBestGeneInPopulation(f), index)
          }
)


test_that("xegaBestinPopulation OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    s<-xegaBestInPopulation(epop10$pop, epop10$fit, lFxegaGaGene)
    expect_identical(s$name==lFxegaGaGene$penv$name(), TRUE)
    expect_equal(s$fitness, max(epop10$fit))
    expect_equal(s$val$fit, s$fit)
    expect_equal(s$phenotype, lFxegaGaGene$DecodeGene(s$genotype, lFxegaGaGene))
          }
)

test_that("xegaBestinPopulation Multiple solutions OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    # we doctor 2 solutions:
    epop10$pop[[1]]$fit<-50
    epop10$pop[[2]]$fit<-50
    cat("epop10$pop[[1]]\n")
    print(epop10$pop[[1]])
    epop10$fit[1]<-50
    epop10$fit[2]<-50
    s<-xegaBestInPopulation(epop10$pop, epop10$fit, lFxegaGaGene, allsolutions=TRUE)
    expect_identical(s$name==lFxegaGaGene$penv$name(), TRUE)
    expect_equal(s$fitness, max(epop10$fit))
    cat("s\n")
    print(s)
    expect_identical((s$value$fit==s$fitness), TRUE)
    expect_equal(s$numberOfSolutions, 2)
    expect_equal(s$phenotype, lFxegaGaGene$DecodeGene(s$genotype, lFxegaGaGene))
    expect_equal(length(s$allgenotypes), 2)
    expect_equal(length(s$allphenotypes), 2)
          }
)

test_that("xegaSummaryPopulation verbose=0  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(0)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    expect_equal(rc, 0)
          }
)

test_that("xegaSummaryPopulation verbose=1  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(1)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 50)
    expect_equal(rc, 0)
          }
)

#
# Snapshot tests?
#

test_that("xegaSummaryPopulation verbose=2  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(2)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 0)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    expect_equal(rc, 0)
          }
)

test_that("xegaSummaryPopulation verbose=3  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(3)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 0)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    expect_equal(rc, 0)
          }
)

test_that("xegaSummaryPopulation verbose=4  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(4)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 0)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    expect_equal(rc, 0)
          }
)

test_that("xegaSummaryPopulation verbose=5  OK",
    {
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    lFxegaGaGene$Verbose<-parm(5)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 0)
    rc<-xegaSummaryPopulation(epop10$pop, epop10$fit, lFxegaGaGene, 5)
    expect_equal(rc, 0)
          }
)

test_that("xegaLogEvalsPopulation  OK",
    {
     evallog<-list()
    pop10<-xegaInitPopulation(10, lFxegaGaGene)
    epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
    evallog<-xegaLogEvalsPopulation(epop10$pop, evallog, 1, lFxegaGaGene)
    evallog<-xegaLogEvalsPopulation(epop10$pop, evallog, 2, lFxegaGaGene)
    expect_equal(evallog[[1]]$fit, epop10$fit[1])    
    expect_equal(evallog[[10]]$fit, epop10$fit[10])    
    expect_equal(evallog[[11]]$fit, epop10$fit[1])    
          }
)

