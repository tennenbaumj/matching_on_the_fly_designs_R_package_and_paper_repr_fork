#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <random>
#include <numeric>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix d_optimal_search_cpp(SEXP P_sexp, int nsim, int n_T) {
    NumericMatrix P_r(P_sexp);
    Eigen::Map<const Eigen::MatrixXd> P(P_r.begin(), P_r.nrow(), P_r.ncol());
    const int n = P.rows();
    const double* p_ptr = P.data();
    const Eigen::VectorXd p_diag = P.diagonal();
    // Global max |P[j,i]| — used for pruning bound: delta(i,j) >= A[i]+B[j] - 2*max_P
    const double max_P = P.cwiseAbs().maxCoeff();
    IntegerMatrix w_mat(n, nsim);

    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    std::random_device rd;
    std::mt19937 g(rd());

    // Hoist per-simulation heap allocations outside the nsim loop.
    Eigen::VectorXd w(n), Pw(n);
    std::vector<int> t_idxs, c_idxs;
    t_idxs.reserve(n_T);
    c_idxs.reserve(n - n_T);

    for (int s = 0; s < nsim; ++s) {
        std::shuffle(indices.begin(), indices.end(), g);
        w.setZero();
        t_idxs.clear();
        c_idxs.clear();
        for (int i = 0; i < n; ++i) {
            if (i < n_T) {
                w(indices[i]) = 1.0;
                t_idxs.push_back(indices[i]);
            } else {
                c_idxs.push_back(indices[i]);
            }
        }

        Pw.noalias() = P * w;
        bool improved = true;
        while (improved) {
            improved = false;
            double best_delta = -1e-10;
            int best_i = -1;
            int best_j = -1;
            int best_t_pos = -1;
            int best_c_pos = -1;

            const double* pw_ptr = Pw.data();

            // Sort t_idxs ascending by A[i] = -2*Pw[i] + p_diag[i]
            // and c_idxs ascending by B[j] = 2*Pw[j] + p_diag[j].
            // With both sorted ascending, A[i]+B[j] is non-decreasing along each axis,
            // enabling inner-j and outer-i early termination via the pruning bound:
            //   delta(i,j) >= A[i]+B[j] - 2*max_P  →  prune when A[i]+B[j] >= best_delta + 2*max_P
            std::sort(t_idxs.begin(), t_idxs.end(), [&](int a, int b) {
                return (-2.0 * pw_ptr[a] + p_diag[a]) < (-2.0 * pw_ptr[b] + p_diag[b]);
            });
            std::sort(c_idxs.begin(), c_idxs.end(), [&](int a, int b) {
                return (2.0 * pw_ptr[a] + p_diag[a]) < (2.0 * pw_ptr[b] + p_diag[b]);
            });

            const double B_min = 2.0 * pw_ptr[c_idxs[0]] + p_diag[c_idxs[0]];

            for (int ti = 0; ti < (int)t_idxs.size(); ++ti) {
                const int i = t_idxs[ti];
                const double Ai = -2.0 * pw_ptr[i] + p_diag[i];
                // Outer-i early termination: if even the smallest B[j] can't unblock this i, stop.
                if (Ai + B_min >= best_delta + 2.0 * max_P) break;

                for (int cj = 0; cj < (int)c_idxs.size(); ++cj) {
                    const int j = c_idxs[cj];
                    const double Bj = 2.0 * pw_ptr[j] + p_diag[j];
                    // Inner-j early termination: A[i]+B[j] is ascending in j; once above threshold, break.
                    if (Ai + Bj >= best_delta + 2.0 * max_P) break;

                    const double delta = Ai + Bj - 2.0 * p_ptr[static_cast<size_t>(j) * n + i];
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_i = i;
                        best_j = j;
                        best_t_pos = ti;
                        best_c_pos = cj;
                        improved = true;
                    }
                }
            }

            if (improved) {
                Pw -= P.col(best_i);
                Pw += P.col(best_j);
                t_idxs[best_t_pos] = best_j;
                c_idxs[best_c_pos] = best_i;
                w(best_i) = 0.0;
                w(best_j) = 1.0;
            }
        }

        for (int i = 0; i < n; ++i) {
            w_mat(i, s) = (int)w(i);
        }
    }

    return w_mat;
}

// [[Rcpp::export]]
IntegerMatrix a_optimal_search_cpp(SEXP P_sexp, SEXP H_sexp, int nsim, int n_T) {
    NumericMatrix P_r(P_sexp);
    NumericMatrix H_r(H_sexp);
    Eigen::Map<const Eigen::MatrixXd> P(P_r.begin(), P_r.nrow(), P_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> H(H_r.begin(), H_r.nrow(), H_r.ncol());
    const int n = P.rows();
    const double* p_ptr = P.data();
    const double* h_ptr = H.data();
    const Eigen::VectorXd p_diag = P.diagonal();
    const Eigen::VectorXd h_diag = H.diagonal();
    // Pruning bound constants: delta_H + obj * delta_P >= C[i]+D[j] - 2*(max_H + obj*max_P)
    const double max_P = P.cwiseAbs().maxCoeff();
    const double max_H = H.cwiseAbs().maxCoeff();
    IntegerMatrix w_mat(n, nsim);

    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    std::random_device rd;
    std::mt19937 g(rd());

    // Hoist per-simulation heap allocations outside the nsim loop.
    Eigen::VectorXd w(n), Pw(n), Hw(n);
    std::vector<int> t_idxs, c_idxs;
    t_idxs.reserve(n_T);
    c_idxs.reserve(n - n_T);

    for (int s = 0; s < nsim; ++s) {
        std::shuffle(indices.begin(), indices.end(), g);
        w.setZero();
        t_idxs.clear();
        c_idxs.clear();
        for (int i = 0; i < n; ++i) {
            if (i < n_T) {
                w(indices[i]) = 1.0;
                t_idxs.push_back(indices[i]);
            } else {
                c_idxs.push_back(indices[i]);
            }
        }

        Pw.noalias() = P * w;
        Hw.noalias() = H * w;
        double wPw = w.dot(Pw);
        double wHw = w.dot(Hw);

        // Objective: (wHw + 1) / (n_T - wPw)
        double obj_curr = (wHw + 1.0) / (n_T - wPw);

        bool improved = true;
        while (improved) {
            improved = false;
            double best_obj = obj_curr - 1e-12;
            int best_i = -1;
            int best_j = -1;
            int best_t_pos = -1;
            int best_c_pos = -1;
            double best_wPw = wPw;
            double best_wHw = wHw;

            const double* pw_ptr = Pw.data();
            const double* hw_ptr = Hw.data();

            // Combined pruning score: C[i] = A_H[i] + obj_curr*A_P[i], D[j] = B_H[j] + obj_curr*B_P[j]
            // Prune pair (i,j) when C[i]+D[j] >= prune_thresh = 2*(max_H + obj_curr*max_P),
            // since then delta_H + obj_curr*delta_P >= 0 for all cross-terms, implying next_obj >= obj_curr.
            // Sort ascending so A[i]+B[j] is non-decreasing along each axis → early termination.
            const double prune_thresh = 2.0 * (max_H + obj_curr * max_P);

            std::sort(t_idxs.begin(), t_idxs.end(), [&](int a, int b) {
                const double Ca = (-2.0*hw_ptr[a] + h_diag[a]) + obj_curr * (-2.0*pw_ptr[a] + p_diag[a]);
                const double Cb = (-2.0*hw_ptr[b] + h_diag[b]) + obj_curr * (-2.0*pw_ptr[b] + p_diag[b]);
                return Ca < Cb;
            });
            std::sort(c_idxs.begin(), c_idxs.end(), [&](int a, int b) {
                const double Da = (2.0*hw_ptr[a] + h_diag[a]) + obj_curr * (2.0*pw_ptr[a] + p_diag[a]);
                const double Db = (2.0*hw_ptr[b] + h_diag[b]) + obj_curr * (2.0*pw_ptr[b] + p_diag[b]);
                return Da < Db;
            });

            // D_min = D[c_idxs[0]] (smallest D, for outer-i early termination)
            const int j0 = c_idxs[0];
            const double D_min = (2.0*hw_ptr[j0] + h_diag[j0]) + obj_curr * (2.0*pw_ptr[j0] + p_diag[j0]);

            for (int ti = 0; ti < (int)t_idxs.size(); ++ti) {
                const int i = t_idxs[ti];
                const double Ci = (-2.0*hw_ptr[i] + h_diag[i]) + obj_curr * (-2.0*pw_ptr[i] + p_diag[i]);
                if (Ci + D_min >= prune_thresh) break;

                for (int cj = 0; cj < (int)c_idxs.size(); ++cj) {
                    const int j = c_idxs[cj];
                    const double Dj = (2.0*hw_ptr[j] + h_diag[j]) + obj_curr * (2.0*pw_ptr[j] + p_diag[j]);
                    if (Ci + Dj >= prune_thresh) break;

                    const double delta_wPw = -2.0 * pw_ptr[i] + 2.0 * pw_ptr[j] +
                        p_diag[i] + p_diag[j] - 2.0 * p_ptr[static_cast<size_t>(j) * n + i];
                    const double delta_wHw = -2.0 * hw_ptr[i] + 2.0 * hw_ptr[j] +
                        h_diag[i] + h_diag[j] - 2.0 * h_ptr[static_cast<size_t>(j) * n + i];

                    const double next_wPw = wPw + delta_wPw;
                    const double next_wHw = wHw + delta_wHw;

                    const double denom = n_T - next_wPw;
                    if (denom <= 1e-10) continue;

                    const double next_obj = (next_wHw + 1.0) / denom;
                    if (next_obj < best_obj) {
                        best_obj = next_obj;
                        best_i = i;
                        best_j = j;
                        best_t_pos = ti;
                        best_c_pos = cj;
                        best_wPw = next_wPw;
                        best_wHw = next_wHw;
                        improved = true;
                    }
                }
            }

            if (improved) {
                obj_curr = best_obj;
                Pw -= P.col(best_i);
                Pw += P.col(best_j);
                Hw -= H.col(best_i);
                Hw += H.col(best_j);
                wPw = best_wPw;
                wHw = best_wHw;
                t_idxs[best_t_pos] = best_j;
                c_idxs[best_c_pos] = best_i;
                w(best_i) = 0.0;
                w(best_j) = 1.0;
            }
        }

        for (int i = 0; i < n; ++i) {
            w_mat(i, s) = (int)w(i);
        }
    }

    return w_mat;
}
