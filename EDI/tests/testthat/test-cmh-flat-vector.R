library(EDI)
library(testthat)

# Reference implementation using R (direct formula, no hash map)
cmh_block_se_ref <- function(y, m_vec, n_total) {
    valid <- !is.na(m_vec) & m_vec > 0 & is.finite(y)
    y <- y[valid]; m_vec <- m_vec[valid]
    blocks <- split(y, m_vec)
    B <- length(blocks)
    if (B == 0) return(NA_real_)
    n_included <- length(y)
    n_b <- n_included / B
    if (n_b <= 1) return(NA_real_)
    var_cmh <- sum(sapply(blocks, function(s_vec) {
        s <- sum(s_vec)
        s * (n_b - s) / (n_b - 1)
    }))
    (2 / n_total) * sqrt(var_cmh)
}

robins_block_se_ref <- function(y, w, m_vec, n_total) {
    valid <- !is.na(m_vec) & m_vec > 0 & is.finite(y) & is.finite(w) & (w == 1 | w == -1)
    y <- y[valid]; w <- w[valid]; m_vec <- m_vec[valid]
    blocks <- split(data.frame(y = y, w = w), m_vec)
    B <- length(blocks)
    if (B == 0 || sum(w == 1) == 0 || sum(w == -1) == 0) return(NA_real_)
    total_t <- sum(w == 1); total_c <- sum(w == -1)
    total_sum_t <- sum(y[w == 1]); total_sum_c <- sum(y[w == -1])
    variance_tot <- sum(sapply(blocks, function(bl) {
        n_b <- nrow(bl)
        if (n_b <= 1) return(NA_real_)
        n_b_half <- n_b / 2
        pH_T <- sum(bl$y[bl$w == 1]) / n_b_half
        pH_C <- sum(bl$y[bl$w == -1]) / n_b_half
        m1 <- max(pH_T, pH_C); m0 <- min(pH_T, pH_C)
        m1*(1-m1)/n_b_half + m0*(1-m0)/n_b_half + ((2*m0-m1)*(1-m1) - m0*(1-m0))/n_b
    })) / B^2
    pH_T <- total_sum_t / total_t; pH_C <- total_sum_c / total_c
    var_rob <- (pH_T*(1-pH_T) + pH_C*(1-pH_C)) / n_total
    sqrt(variance_tot + var_rob)
}

make_cmh_data <- function(n, B, seed = 42) {
    set.seed(seed)
    y <- rbinom(n, 1, 0.4)
    m_vec <- rep(seq_len(B), each = n / B)
    list(y = as.numeric(y), m_vec = as.integer(m_vec))
}

make_robins_data <- function(n, B, seed = 42) {
    set.seed(seed)
    y <- rbinom(n, 1, 0.4)
    m_vec <- rep(seq_len(B), each = n / B)
    w <- rep(c(1, -1), n / 2)
    list(y = as.numeric(y), w = as.numeric(w), m_vec = as.integer(m_vec))
}

# --- compute_cmh_block_se_cpp tests ---

test_that("cmh: basic result matches R reference", {
    d <- make_cmh_data(200, 50)
    cpp_val <- EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, 200L)
    ref_val <- cmh_block_se_ref(d$y, d$m_vec, 200)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("cmh: small case (n=20, B=5) matches reference", {
    d <- make_cmh_data(20, 5)
    cpp_val <- EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, 20L)
    ref_val <- cmh_block_se_ref(d$y, d$m_vec, 20)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("cmh: larger case (n=1000, B=200) matches reference", {
    d <- make_cmh_data(1000, 200)
    cpp_val <- EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, 1000L)
    ref_val <- cmh_block_se_ref(d$y, d$m_vec, 1000)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("cmh: all y=0 in some blocks handled correctly", {
    y <- c(0, 0, 1, 0, 1, 0, 0, 0)
    m_vec <- as.integer(c(1, 1, 2, 2, 3, 3, 4, 4))
    cpp_val <- EDI:::compute_cmh_block_se_cpp(as.numeric(y), m_vec, 8L)
    ref_val <- cmh_block_se_ref(as.numeric(y), m_vec, 8)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("cmh: non-contiguous block IDs handled correctly", {
    y <- c(1, 0, 1, 1, 0, 1)
    m_vec <- as.integer(c(1, 1, 3, 3, 7, 7))  # gaps in IDs
    cpp_val <- EDI:::compute_cmh_block_se_cpp(as.numeric(y), m_vec, 6L)
    ref_val <- cmh_block_se_ref(as.numeric(y), m_vec, 6)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("cmh: returns NA when n_total <= 0", {
    d <- make_cmh_data(20, 5)
    expect_true(is.na(EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, 0L)))
    expect_true(is.na(EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, -1L)))
})

test_that("cmh: returns NA for all-NA m_vec", {
    y <- c(1, 0, 1, 0)
    m_vec <- as.integer(c(NA, NA, NA, NA))
    result <- EDI:::compute_cmh_block_se_cpp(as.numeric(y), m_vec, 4L)
    expect_true(is.na(result))
})

test_that("cmh: result is positive", {
    d <- make_cmh_data(100, 25, seed = 7)
    cpp_val <- EDI:::compute_cmh_block_se_cpp(d$y, d$m_vec, 100L)
    expect_true(is.finite(cpp_val) && cpp_val > 0)
})

# --- compute_extended_robins_block_se_cpp tests ---

test_that("robins: basic result matches R reference", {
    d <- make_robins_data(200, 50)
    cpp_val <- EDI:::compute_extended_robins_block_se_cpp(d$y, d$w, d$m_vec, 200L)
    ref_val <- robins_block_se_ref(d$y, d$w, d$m_vec, 200)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("robins: small case (n=20, B=5) matches reference", {
    d <- make_robins_data(20, 5)
    cpp_val <- EDI:::compute_extended_robins_block_se_cpp(d$y, d$w, d$m_vec, 20L)
    ref_val <- robins_block_se_ref(d$y, d$w, d$m_vec, 20)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("robins: larger case (n=1000, B=200) matches reference", {
    d <- make_robins_data(1000, 200)
    cpp_val <- EDI:::compute_extended_robins_block_se_cpp(d$y, d$w, d$m_vec, 1000L)
    ref_val <- robins_block_se_ref(d$y, d$w, d$m_vec, 1000)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("robins: non-contiguous block IDs handled correctly", {
    set.seed(5)
    y <- rbinom(12, 1, 0.4)
    m_vec <- as.integer(rep(c(1, 3, 5, 10), each = 3))
    w <- as.numeric(rep(c(1, -1, 1), 4))
    cpp_val <- EDI:::compute_extended_robins_block_se_cpp(as.numeric(y), w, m_vec, 12L)
    ref_val <- robins_block_se_ref(as.numeric(y), w, m_vec, 12)
    expect_equal(cpp_val, ref_val, tolerance = 1e-10)
})

test_that("robins: returns NA when n_total <= 0", {
    d <- make_robins_data(20, 5)
    expect_true(is.na(EDI:::compute_extended_robins_block_se_cpp(d$y, d$w, d$m_vec, 0L)))
})

test_that("robins: result is positive", {
    d <- make_robins_data(100, 25, seed = 9)
    cpp_val <- EDI:::compute_extended_robins_block_se_cpp(d$y, d$w, d$m_vec, 100L)
    expect_true(is.finite(cpp_val) && cpp_val > 0)
})
