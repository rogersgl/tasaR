
test_that("batchAnalyzeAmplicon accepts different input types and returns the expected results", {
  expect_output(amp.ctrl <- analyzeAmplicon(test.settings))
  expect_output(amp.batch.r <- batchAnalyzeAmplicon(rbind(make.set.df(), make.set.df())))
  expect_output(amp.batch.f <- batchAnalyzeAmplicon(file.path(tempdir(), "settings-batch.csv")))

  expect_equal(class(amp.batch.r), "list")
  expect_equal(class(amp.batch.f), "list")

  for (i in seq_along(amp.batch.r)) {
    expect_equal(amp.batch.r[[i]], amp.ctrl)
  }
  for (i in seq_along(amp.batch.f)) {
    expect_equal(amp.batch.f[[i]], amp.ctrl)
  }
})


test_that("batchSummarize obeys flags correctly and returns the expected data", {
  expect_output(results.list <- batchAnalyzeAmplicon(rbind(make.set.df(), make.set.df())))
  expect_visible(batchSummarize(results.list))
  expect_invisible(batchSummarize(results.list, suppressConsoleOutput = TRUE))
  expect_warning(batchSummarize(results.list, export = TRUE, path = NULL))
  suppressMessages(test.summary.list <- batchSummarize(results.list, export = TRUE, path = file.path(tempdir(), "summary")))
  expect_equal(class(test.summary.list), "list")
  expect_equal(length(test.summary.list), 2)
  expect_equal(names(test.summary.list), c("Graphs", "Tables"))
  expect_equal(names(test.summary.list$Graphs), c("AvgMutAll", "AvgMutCyt", "MutPosDNA", "MutPosAA", "MutCytHM", "AIDBox", "MutTypes"))
  expect_all_true(sapply(test.summary.list$Graphs, function(x) {
    "ggplot" %in% class(x)
  }))
  expect_equal(names(test.summary.list$Tables), c("MotifSums", "AllMutationsDNA", "CytosineMutations", "NonCytosineMutations", "AllMutationsAA", "MutationTypes"))
  expect_true(file.exists(file.path(tempdir(), "summary", "Motif Sums (batch).csv")))
  expect_true(file.exists(file.path(tempdir(), "summary", "All DNA Mutations (batch).csv")))
  expect_true(file.exists(file.path(tempdir(), "summary", "Cytosine DNA Mutations (batch).csv")))
  expect_true(file.exists(file.path(tempdir(), "summary", "Non-Cytosine DNA Mutations (batch).csv")))
  expect_true(file.exists(file.path(tempdir(), "summary", "All Protein Mutations (batch).csv")))
  expect_true(file.exists(file.path(tempdir(), "summary", "DNA Repair Types (batch).csv")))

  expect_true(file.exists(file.path(tempdir(), "summary", "Average DNA Mutation Rate.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "Average AID Cytosine Mutation Rate.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "DNA Mutation Distributions.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "Protein Mutation Distributions.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "AID Cytosine Mutation Heatmap.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "AID Mutation Boxplot.pdf")))
  expect_true(file.exists(file.path(tempdir(), "summary", "DNA Repair Types.pdf")))
})
