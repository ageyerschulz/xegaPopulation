library(testthat)
library(xegaPopulation)

test_that("xegaConfiguration OK",
          { BG<-7
            XORBNF<-8
  GA<-function(pe, gr=NULL, nope=1.5, sle="test", ok=TRUE)
  {xegaConfiguration("GA", substitute(pe), substitute(gr), environment())}
	  a<-GA(BG)
	   expect_setequal(names(a), c("GAconf", "GAenv"))
	   expect_identical(a$GAenv$pe, 7)
	   expect_identical(a$GAenv$gr, NULL)
	   expect_identical(a$GAenv$nope, 1.5)
	   expect_identical(a$GAenv$sle, "test")
	   expect_identical(a$GAenv$ok, TRUE)
            WHO<-10
          b<-GA(pe=WHO, gr=XORBNF, nope=0.0, sle=GA, ok=FALSE)
	   expect_identical(b$GAenv$pe, 10)
	   expect_identical(b$GAenv$gr, 8)
	   expect_identical(b$GAenv$nope, 0.0)
	   expect_identical(b$GAenv$sle, GA)
	   expect_identical(b$GAenv$ok, FALSE)
          }
)


