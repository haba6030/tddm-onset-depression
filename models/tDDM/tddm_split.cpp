// Split-at-tau trial likelihood in C++ (2026-09-28). Same algorithm as the R engine's
// accurate path in tddm_engine.R: QUADPACK Rdqags (the routine behind R integrate()) with
// rel.tol = abs.tol = 1e-6, 1000 subdivisions, breakpoints {0, h, bsep-h, bsep} with
// h = min(bsep/4, 20 sqrt(rem)), and the same 4001-point trapezoid fallback on any error.
// Verified against the R path: scripts/tDDM/verify_cpp_likelihood.R.
#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <cmath>
using namespace Rcpp;
static const int NIMAGE = 4;
static inline double dw_small(double t, double v, double a, double w) {
  double pre = (a / std::sqrt(2 * M_PI * t * t * t)) * std::exp(-v * a * w - v * v * t / 2), s = 0;
  for (int k = -7; k <= 7; k++) { double wk = w + 2 * k; s += wk * std::exp(-wk * wk * a * a / (2 * t)); }
  return pre * s;
}
static inline double dw_large(double t, double v, double a, double w) {
  double pre = (M_PI / (a * a)) * std::exp(-v * a * w - v * v * t / 2), s = 0;
  for (int k = 1; k <= 7; k++) s += k * std::sin(k * M_PI * w) * std::exp(-k * k * M_PI * M_PI * t / (2 * a * a));
  return pre * s;
}
static inline double dwiener(double t, double v, double a, double w) {
  if (t <= 0 || a <= 0) return 0;
  return (t < a * a / (M_PI * M_PI)) ? dw_small(t, v, a, w) : dw_large(t, v, a, w);
}
static inline double trans_image(double x, double t, double x0, double v, double a) {
  if (t <= 0) return 0;
  double s = 0;
  for (int n = -NIMAGE; n <= NIMAGE; n++)
    s += std::exp(-(x - x0 - 2 * n * a) * (x - x0 - 2 * n * a) / (2 * t)) - std::exp(-(x + x0 - 2 * n * a) * (x + x0 - 2 * n * a) / (2 * t));
  if (s <= 0) return 0;
  double val = std::exp(v * (x - x0) - v * v * t / 2 - 0.5 * std::log(2 * M_PI * t)) * s;
  return (std::isfinite(val) && val > 0) ? val : 0;
}
struct Pars { double d1, d2, tau, bsep, w, rem; int ch; bool bad; };
static void integrand(double *x, int n, void *ex) {
  Pars *p = (Pars *) ex;
  for (int i = 0; i < n; i++) {
    double pr = trans_image(x[i], p->tau, p->bsep * p->w, p->d1, p->bsep);
    double pf = (p->ch == 1) ? dwiener(p->rem, -p->d2, p->bsep, 1 - x[i] / p->bsep) : dwiener(p->rem, p->d2, p->bsep, x[i] / p->bsep);
    x[i] = pr * pf;
    if (!std::isfinite(x[i])) p->bad = true;   // R integrate() stops on non-finite values
  }
}
// [[Rcpp::export]]
double ll_trial_cpp(double d1, double d2, double tau, double bsep, double w, double nDT, double rt, int ch) {
  double dt = rt - nDT;
  if (dt <= 0.001) return 1e-12;
  if (dt < tau || tau <= 0) {
    double eff = (tau <= 0) ? d2 : d1;
    double d = (ch == 1) ? dwiener(dt, -eff, bsep, 1 - w) : dwiener(dt, eff, bsep, w);
    return std::max(d, 1e-12);
  }
  Pars p = {d1, d2, tau, bsep, w, dt - tau, ch, false};
  double h = std::min(bsep / 4, 20 * std::sqrt(p.rem));
  double br[4] = {0, h, bsep - h, bsep}; double tot = 0; bool fail = false;
  for (int k = 0; k < 3 && !fail; k++) {
    if (br[k + 1] <= br[k]) continue;                // R: unique() drops repeated breakpoints
    double a = br[k], b = br[k + 1], epsabs = 1e-6, epsrel = 1e-6, res = 0, abserr = 0;
    int neval = 0, ier = 0, limit = 1000, lenw = 4 * limit, last = 0;
    std::vector<int> iwork(limit); std::vector<double> work(lenw);
    Rdqags(integrand, &p, &a, &b, &epsabs, &epsrel, &res, &abserr, &neval, &ier, &limit, &lenw, &last, iwork.data(), work.data());
    if (ier != 0 || p.bad) fail = true;              // R integrate() errors on any ier > 0
    tot += res;
  }
  if (fail) {                                        // same fallback as the R engine
    int n = 4001; double dx = bsep / (n - 1); std::vector<double> xg(n);
    for (int i = 0; i < n; i++) xg[i] = i * dx;
    p.bad = false; integrand(xg.data(), n, &p);
    tot = 0; for (int i = 0; i < n - 1; i++) tot += (xg[i] + xg[i + 1]) / 2;
    tot *= dx;
  }
  return std::max(tot, 1e-12);
}
// [[Rcpp::export]]
double nll_tau_cpp(NumericVector par, NumericVector s1, NumericVector s2, NumericVector rt, IntegerVector ch, double sigma, double eps) {
  double bsep = 2 * par[2] / sigma; if (bsep <= 0.01 || par[3] <= 0) return 1e12;
  double nll = 0;
  for (int i = 0; i < rt.size(); i++) {
    double d1 = par[0] * s1[i] / sigma, d2 = (par[0] * s1[i] + par[1] * s2[i]) / sigma;
    double l = ll_trial_cpp(d1, d2, par[4], bsep, 0.5, par[3], rt[i], ch[i]);
    if (ISNAN(l)) return NA_REAL;
    nll -= std::log(std::max(l, eps));
  }
  return std::isfinite(nll) ? nll : 1e12;
}
