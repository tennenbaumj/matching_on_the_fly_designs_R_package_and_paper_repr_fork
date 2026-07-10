#include <RcppEigen.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Internal result of evaluating p and dp/d(delta) at a single point.
struct LrtEval {
	double p;
	double dp;
	bool valid;
	LrtEval() : p(R_NaReal), dp(R_NaReal), valid(false) {}
};

// Evaluate the LRT p-value and its derivative at delta via R callbacks.
// Envelope theorem: dp/d(delta) = 2 * f_chi2(T) * score[j-1]
static LrtEval eval_lrt(
	Function& fit_null_fn,
	Function& neg_loglik_fn,
	Function& score_fn,
	double delta,
	double full_negloglik,
	int j   // 1-indexed treatment position in score vector
) {
	LrtEval res;

	SEXP null_fit_sexp;
	try {
		null_fit_sexp = fit_null_fn(delta);
	} catch (...) {
		return res;
	}
	if (Rf_isNull(null_fit_sexp)) return res;
	RObject null_fit(null_fit_sexp);

	double null_nll;
	try {
		SEXP nll_sexp = neg_loglik_fn(null_fit);
		if (Rf_isNull(nll_sexp)) return res;
		null_nll = as<double>(nll_sexp);
	} catch (...) {
		return res;
	}
	if (!R_finite(null_nll)) return res;

	double T_stat = std::max(2.0 * (null_nll - full_negloglik), 0.0);
	res.p = R::pchisq(T_stat, 1.0, 0, 0);   // lower.tail = FALSE
	res.valid = true;

	// Derivative via envelope theorem (optional — missing dp just forces bisection)
	try {
		SEXP sc_sexp = score_fn(null_fit);
		if (!Rf_isNull(sc_sexp)) {
			NumericVector sc(sc_sexp);
			if (j >= 1 && j <= (int)sc.size()) {
				double score_j = sc[j - 1];
				if (T_stat > 0.0 && R_finite(score_j)) {
					res.dp = 2.0 * R::dchisq(T_stat, 1.0, 0) * score_j;
				}
			}
		}
	} catch (...) {
		// dp stays NA; loop falls back to pure bisection step
	}

	return res;
}

//' LRT confidence interval by Newton-Raphson + bisection (Rcpp implementation)
//'
//' Implements the bracket search + NR+bisection loop entirely in C++, calling
//' back into R only for \code{fit_null_fn}, \code{neg_loglik_fn}, and
//' \code{score_fn}. The derivative \eqn{dp/d\delta = 2 f_{\chi^2}(T) \cdot
//' \mathrm{score}[j]} follows from the envelope theorem.
//'
//' @param fit_null_fn   R function \code{delta -> list} (constrained null fit)
//' @param neg_loglik_fn R function \code{fit -> double} (negative log-likelihood)
//' @param score_fn      R function \code{fit -> numeric} (score vector)
//' @param est           Point estimate of the treatment effect
//' @param full_negloglik Negative log-likelihood of the unrestricted model
//' @param alpha         Significance level (e.g. 0.05)
//' @param step          Initial step size for exponential bracket search
//' @param lower_seed    Initial lower-bound candidate, typically the Wald lower CI
//' @param upper_seed    Initial upper-bound candidate, typically the Wald upper CI
//' @param j             1-indexed position of the treatment coefficient in the
//'   score vector
//' @param max_bracket   Maximum exponential bracket search iterations (default 60)
//' @param max_nr_iter   Maximum NR+bisection iterations per bound (default 25)
//' @param tol_p         P-value convergence tolerance (default 1e-7)
//' @param tol_bracket   Bracket-width convergence tolerance (default 1e-8)
//'
//' @return Unnamed numeric vector of length 2: \code{[lower_bound, upper_bound]}
//'
// [[Rcpp::export]]
NumericVector lrt_ci_nr_cpp(
	Function fit_null_fn,
	Function neg_loglik_fn,
	Function score_fn,
	double est,
	double full_negloglik,
	double alpha,
	double step,
	double lower_seed,
	double upper_seed,
	int j,
	int max_bracket = 60,
	int max_nr_iter = 25,
	double tol_p     = 1e-7,
	double tol_bracket = 1e-8
) {
	NumericVector ci(2, R_NaReal);

	for (int dir_idx = 0; dir_idx < 2; ++dir_idx) {
		double direction = (dir_idx == 0) ? -1.0 : 1.0;
		double seed = (dir_idx == 0) ? lower_seed : upper_seed;

		// Phase 1: begin at the Wald-seeded candidate, then expand outward if needed.
		double delta_outer = R_NaReal;
		LrtEval ev_outer;
		if (R_finite(seed) && ((direction < 0.0 && seed < est) || (direction > 0.0 && seed > est))) {
			LrtEval ev_seed = eval_lrt(fit_null_fn, neg_loglik_fn, score_fn, seed, full_negloglik, j);
			if (ev_seed.valid && R_finite(ev_seed.p) && ev_seed.p <= alpha) {
				delta_outer = seed;
				ev_outer = ev_seed;
			}
		}
		if (!R_finite(delta_outer)) {
			for (int i = 0; i < max_bracket; ++i) {
				double d = est + direction * step * std::pow(2.0, (double)i);
				LrtEval ev = eval_lrt(fit_null_fn, neg_loglik_fn, score_fn, d, full_negloglik, j);
				if (ev.valid && R_finite(ev.p) && ev.p <= alpha) {
					delta_outer = d;
					ev_outer = ev;
					break;
				}
			}
		}
		if (!R_finite(delta_outer)) continue;

		// Phase 2: NR + bisection
		// Invariant: p(a) > alpha, p(b) <= alpha
		double a = est;
		double b = delta_outer;
		double delta = b;
		LrtEval ev = ev_outer;

		bool converged_p = false;
		for (int k = 0; k < max_nr_iter; ++k) {
			if (ev.valid && std::abs(ev.p - alpha) < tol_p) {
				converged_p = true;
				break;
			}
			double lo = std::min(a, b);
			double hi = std::max(a, b);
			if (hi - lo < tol_bracket) break;

			// Newton step; fall back to bisection if dp is zero/NA or step exits bracket
			double delta_nr = R_NaReal;
			if (ev.valid && R_finite(ev.dp) && std::abs(ev.dp) > 1e-15) {
				delta_nr = delta - (ev.p - alpha) / ev.dp;
			}

			double delta_new;
			if (R_finite(delta_nr) && delta_nr > lo && delta_nr < hi) {
				delta_new = delta_nr;
			} else {
				delta_new = (lo + hi) / 2.0;
			}

			LrtEval ev_new = eval_lrt(fit_null_fn, neg_loglik_fn, score_fn, delta_new, full_negloglik, j);
			if (!ev_new.valid || !R_finite(ev_new.p)) {
				delta_new = (lo + hi) / 2.0;
				ev_new = eval_lrt(fit_null_fn, neg_loglik_fn, score_fn, delta_new, full_negloglik, j);
				if (!ev_new.valid || !R_finite(ev_new.p)) break;
			}

			if (ev_new.p > alpha) a = delta_new; else b = delta_new;
			delta = delta_new;
			ev = ev_new;
		}

		ci[dir_idx] = converged_p ? delta : (a + b) / 2.0;
	}

	return ci;
}

//' CI inversion by p-value bracket search + bisection (Rcpp implementation)
//'
//' Used by \code{invert_test_pval_confidence_interval} (score CI and any other
//' CI that inverts a scalar p-value function with no derivative available).
//' The Wald bounds are tried first as bracket candidates before falling back
//' to exponential search.  Bisection then polishes to \code{tol} in delta-space.
//'
//' @param pval_fn     R function \code{delta -> double} two-sided p-value
//' @param est         Point estimate of the treatment effect
//' @param alpha       Significance level (e.g. 0.05)
//' @param step        Initial step size for exponential bracket search
//' @param lower_seed  Wald lower CI bound (used as first bracket candidate;
//'   pass \code{NA_real_} to skip)
//' @param upper_seed  Wald upper CI bound (same)
//' @param max_bracket Maximum exponential bracket search iterations (default 60)
//' @param max_bisect  Maximum bisection iterations (default 60)
//' @param tol         Bracket-width convergence tolerance in delta-space (default 1e-6)
//'
//' @return Unnamed numeric vector of length 2: \code{[lower_bound, upper_bound]}
//'
// [[Rcpp::export]]
NumericVector pval_invert_ci_cpp(
	Function pval_fn,
	double est,
	double alpha,
	double step,
	double lower_seed,
	double upper_seed,
	int max_bracket = 60,
	int max_bisect  = 60,
	double tol      = 1e-6
) {
	NumericVector ci(2, R_NaReal);

	// Evaluate p at the estimate once — used to guard degenerate cases
	double p_est = R_NaReal;
	try {
		SEXP s = pval_fn(est);
		if (!Rf_isNull(s)) p_est = as<double>(s);
	} catch (...) {}
	if (!R_finite(p_est)) return ci;
	// If p at the estimate is already below alpha the CI degenerates to a point
	if (p_est < alpha) { ci[0] = est; ci[1] = est; return ci; }

	auto eval_p = [&](double delta) -> double {
		try {
			SEXP s = pval_fn(delta);
			if (Rf_isNull(s)) return R_NaReal;
			return as<double>(s);
		} catch (...) {
			return R_NaReal;
		}
	};

	for (int dir_idx = 0; dir_idx < 2; ++dir_idx) {
		double direction = (dir_idx == 0) ? -1.0 : 1.0;
		double seed = (dir_idx == 0) ? lower_seed : upper_seed;

		// Phase 1: bracket search
		// Invariant we need: p(a) > alpha (inner), p(b) <= alpha (outer)
		double a = est;   // always inner (p >= alpha)
		double b = R_NaReal;

		// Try Wald seed first
		if (R_finite(seed) && ((direction < 0.0 && seed < est) || (direction > 0.0 && seed > est))) {
			double p_seed = eval_p(seed);
			if (R_finite(p_seed) && p_seed <= alpha) {
				b = seed;
			}
		}

		// Exponential fallback
		if (!R_finite(b)) {
			for (int i = 0; i < max_bracket; ++i) {
				double d = est + direction * step * std::pow(2.0, (double)i);
				double p_d = eval_p(d);
				if (R_finite(p_d) && p_d <= alpha) {
					b = d;
					break;
				}
			}
		}
		if (!R_finite(b)) continue;

		// Phase 2: bisection
		// [a, b] bracket: p(a) > alpha, p(b) <= alpha
		for (int k = 0; k < max_bisect; ++k) {
			double lo = std::min(a, b);
			double hi = std::max(a, b);
			if (hi - lo < tol) break;
			double mid = (lo + hi) / 2.0;
			double p_mid = eval_p(mid);
			if (!R_finite(p_mid)) break;
			if (p_mid > alpha) a = mid; else b = mid;
		}

		ci[dir_idx] = b;
	}

	return ci;
}
