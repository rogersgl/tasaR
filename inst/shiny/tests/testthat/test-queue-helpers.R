test_that("sample order reconciliation preserves existing order and appends new ids", {
  expect_equal(
    .reconcile_sample_order(c("a", "b", "c"), c("b", "z", "a", "b")),
    c("b", "a", "c")
  )
})

test_that("sample data frames can be ordered and moved", {
  df <- data.frame(
    id = c("a", "b", "c"),
    value = 1:3,
    check.names = FALSE
  )

  expect_equal(.order_sample_df(df, c("c", "a"), "id")$id, c("c", "a", "b"))
  expect_equal(.move_order_item(c("a", "b", "c"), "b", "up"), c("b", "a", "c"))
  expect_equal(.move_order_item(c("a", "b", "c"), "b", "down"), c("a", "c", "b"))
  expect_equal(.move_order_item(c("a", "b", "c"), "a", "up"), c("a", "b", "c"))
})

test_that("rows can be removed by id", {
  df <- data.frame(
    "Queue ID" = c("a", "b"),
    "Sample Name" = c("sample_a", "sample_b"),
    check.names = FALSE
  )

  res <- .remove_batch_row_by_id(df, "b", id_col = "Queue ID")
  expect_true(res$success)
  expect_equal(res$removed[["Sample Name"]], "sample_b")
  expect_equal(res$df[["Queue ID"]], "a")

  missing <- .remove_batch_row_by_id(df, "z", id_col = "Queue ID")
  expect_false(missing$success)
  expect_equal(missing$df, df)
})

test_that("queue action HTML keeps expected classes and labels", {
  html <- .queue_actions_html("row_1", "remove_id", "move_up_id", "move_down_id")

  expect_match(html, "queue-move-btn", fixed = TRUE)
  expect_match(html, "queue-remove-btn", fixed = TRUE)
  expect_match(html, "aria-label='Move up'", fixed = TRUE)
  expect_match(html, "aria-label='Move down'", fixed = TRUE)
  expect_match(html, "aria-label='Remove'", fixed = TRUE)
})
