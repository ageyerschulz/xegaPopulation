library(testthat)
library(xegaSelectGene)
library(xegaGaGene)
library(xegaPopulation)

test_that("ExponentialMultiplicativeCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
                   Alpha=parm(0.99))
          expect_identical(ExponentialMultiplicativeCooling(0, lF), 100)
          expect_equal(ExponentialMultiplicativeCooling(1, lF), 99)
          expect_equal(ExponentialMultiplicativeCooling(5, lF), (100*0.99^5))
          }
)

test_that("ExponentialMultiplicativeCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
                   Alpha=parm(0.99))
          expect_identical(LogarithmicMultiplicativeCooling(0, lF), 100)
          expect_equal(LogarithmicMultiplicativeCooling(1, lF), 59.30439,
	  tolerance=0.0001)
          expect_equal(LogarithmicMultiplicativeCooling(5, lF), 36.05108,
	  tolerance=0.0001)
          }
)

test_that("PowerMultiplicativeCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
                   Alpha=parm(0.99),
	           CoolingPower=parm(2))
          expect_identical(PowerMultiplicativeCooling(0, lF), 100)
          expect_equal(PowerMultiplicativeCooling(1, lF), 50.25126,
	  tolerance=0.0001)
          expect_equal(PowerMultiplicativeCooling(5, lF), 3.883495,
	  tolerance=0.0001)
          }
)

test_that("PowerAdditiveCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
		   TempN=parm(10),
                   Generations=parm(50),
	           CoolingPower=parm(2))
          expect_identical(PowerAdditiveCooling(0, lF), 100)
          expect_equal(PowerAdditiveCooling(1, lF), 96.436,
	  tolerance=0.0001)
          expect_equal(PowerAdditiveCooling(2, lF), 92.94,
	  tolerance=0.0001)
          expect_equal(PowerAdditiveCooling(50, lF), 10,
	  tolerance=0.0001)
          }
)

test_that("ExponenentialAdditiveCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
		   TempN=parm(10),
                   Generations=parm(50))
          expect_equal(ExponentialAdditiveCooling(0, lF), 99.01099,
	  tolerance=0.0001)
          expect_equal(ExponentialAdditiveCooling(1, lF), 98.81851,
	  tolerance=0.0001)
          expect_equal(ExponentialAdditiveCooling(2, lF), 98.58916,
	  tolerance=0.0001)
          expect_equal(ExponentialAdditiveCooling(50, lF), 10.98901,
	  tolerance=0.0001)
          }
)

test_that("TrigonometricAdditiveCooling OK",
          {
          parm<-function(x){function() {return(x)}}
          lF<-list(Temp0=parm(100),
		   TempN=parm(10),
                   Generations=parm(50))
          expect_equal(TrigonometricAdditiveCooling(0, lF), 100,
	  tolerance=0.0001)
          expect_equal(TrigonometricAdditiveCooling(1, lF), 99.9112,
	  tolerance=0.0001)
          expect_equal(TrigonometricAdditiveCooling(2, lF), 99.64516,
	  tolerance=0.0001)
          expect_equal(TrigonometricAdditiveCooling(50, lF), 10,
	  tolerance=0.0001)
          }
)

test_that("CoolingFactory OK",
          {
          Fun<-CoolingFactory()
          expect_identical(body(Fun), body(ExponentialMultiplicativeCooling))
          Fun<-CoolingFactory("ExponentialMultiplicative")
          expect_identical(body(Fun), body(ExponentialMultiplicativeCooling))
          Fun<-CoolingFactory("LogarithmicMultiplicative")
          expect_identical(body(Fun), body(LogarithmicMultiplicativeCooling))
          Fun<-CoolingFactory("PowerMultiplicative")
          expect_identical(body(Fun), body(PowerMultiplicativeCooling))
          Fun<-CoolingFactory("PowerAdditive")
          expect_identical(body(Fun), body(PowerAdditiveCooling))
          Fun<-CoolingFactory("ExponentialAdditive")
          expect_identical(body(Fun), body(ExponentialAdditiveCooling))
          Fun<-CoolingFactory("TrigonometricAdditive")
          expect_identical(body(Fun), body(TrigonometricAdditiveCooling))
          expect_error(CoolingFactory("Stchastic"))
          }
)

