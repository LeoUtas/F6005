#include <TMB.hpp>
// #include <math.h>

template <class Type>
Type power(Type x, int y) // x is defined as Type for later use of unspecified x (line 34)
{
  Type result = 1;
  for (int i = 0; i < y; i++)
  {
    result = result * x;
  }
  return result;
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(ssb);
  DATA_VECTOR(rec);
  DATA_VECTOR(log_ssb);
  DATA_VECTOR(log_rec);

  PARAMETER(log_alpha);
  PARAMETER(log_beta);
  PARAMETER(log_sd_log_rec_me);
  PARAMETER(logit_ar_log_rec_me);

  int n = ssb.size();

  Type alpha = exp(log_alpha);
  Type beta = exp(log_beta);
  Type sd_log_rec_me = exp(log_sd_log_rec_me);
  Type ar_log_rec_me = exp(logit_ar_log_rec_me) / (1.0 + exp(logit_ar_log_rec_me));
  Type st_sd_log_rec_me = sd_log_rec_me / sqrt(1.0 - power(ar_log_rec_me, 2));

  Type nll = 0.0;
  Type zero = 0.0;

  vector<Type> rec_pred = alpha * ssb / (beta + ssb);
  vector<Type> log_rec_pred = log(rec_pred);

  vector<Type> log_rec_resid = log_rec - log_rec_pred;
  vector<Type> log_rec_resid_std = log_rec_resid / st_sd_log_rec_me;

  // nll for log recruitment;
  nll -= dnorm(log_rec_resid(0), zero, st_sd_log_rec_me, true);
  // L34 does not work if 0.0 is used instead of zero??;
  for (int i = 1; i < n; ++i)
  {
    nll -= dnorm(log_rec_resid(i), ar_log_rec_me * log_rec_resid(i - 1), sd_log_rec_me, true);
  }

  //  using namespace density;
  //  nll += SCALE(AR1(ar_lrec_me),st_sd_lrec_me)(lrec_resid);

  REPORT(alpha);
  REPORT(beta);
  REPORT(sd_log_rec_me);
  REPORT(ar_log_rec_me);
  REPORT(st_sd_log_rec_me);

  REPORT(rec_pred);
  REPORT(log_rec_pred);
  REPORT(log_rec_resid);
  REPORT(log_rec_resid_std);

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sd_log_rec_me);
  ADREPORT(ar_log_rec_me);
  ADREPORT(st_sd_log_rec_me);

  return nll;
}
