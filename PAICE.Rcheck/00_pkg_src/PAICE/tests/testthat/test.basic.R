library(PAICE)
context("Basic test")

test_that("Testing colonization inference", {
     expect_equal({
          x <- colonization(CmonsData, CmonsNetwork)
          x$Total[1, 1]
          }, 13)
})

test_that("Testing maximum of colonization events with a few replicates", {
     expect_equal({
          set.seed(24)
          rare <- rarecol(data               = CmonsData,
                          network            = CmonsNetwork,
                          replicates_field   = 5,
                          replicates_genetic = 5,
                          monitor            = FALSE)
          max <- maxCol(rare)
          c(round(x      = max$Summary[1, 1],
                  digits = 4),
            round(x      = max$Summary[2, 1],
                  digits = 4))
     }, c(34.0653, 41.1658))
})
