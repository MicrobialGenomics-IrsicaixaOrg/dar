# Function bake() @ L85

    Code
      testthat::expect_s3_class(cool(test), "tbl_df")
    Message
      i Baking with count_cutoff = 1
      Excluding "maaslin__Gundain".

---

    Code
      testthat::expect_equal(nrow(cool(test)), 27)
    Message
      i Baking with count_cutoff = 1
      Excluding "maaslin__Gundain".

---

    Code
      testthat::expect_s3_class(cool(test), "tbl_df")
    Message
      i Baking with count_cutoff = 2

---

    Code
      testthat::expect_equal(nrow(cool(test)), 41)
    Message
      i Baking with count_cutoff = 2

