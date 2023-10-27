library(PAICE)
context("NAMESPACE test")

test_that("change basic functions", {
     expect_equal({
          aggregate <- function(x) 2
          x <- colonization(CmonsData, CmonsNetwork)
          x$Total[1, 1]
     }, 13)
})

test_that("change package functions", {
     expect_equal({
          colonization <- function(x) x
          x1 <- colonization(2)
          x2 <- PAICE::colonization(CmonsData, CmonsNetwork)
          x <- c(x1, x2$Total[1, 1])
          rm(colonization)
          y <- colonization(CmonsData, CmonsNetwork)
          x <- c(x, y$Total[1, 1])
          x
     }, c(2, 13, 13))
     expect_equal({
          x <- rarecol(data               = CmonsData,
                       network            = CmonsNetwork,
                       replicates_field   = 1,
                       replicates_genetic = 1,
                       monitor            = FALSE)
          fun <- function(x) stop("error")
          x <- maxCol(x)
          rm(fun)
          class(x)
     }, "maxCol")
})
