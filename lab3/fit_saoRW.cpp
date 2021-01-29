#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(ssb);
  DATA_VECTOR(rec);
  DATA_VECTOR(log_ssb);
  DATA_VECTOR(log_rec);

  PARAMETER(log_sao);
  PARAMETER(log_beta);
  PARAMETER(log_sd_log_rec_me);
  PARAMETER(log_sd_log_sao_dev);

  PARAMETER_VECTOR(log_sao_dev);

  int n = ssb.size();

  Type beta = exp(log_beta);
  vector<Type> alpha = beta * exp(log_sao + log_sao_dev);
  Type sd_log_rec_me = exp(log_sd_log_rec_me);
  Type sd_log_sao_dev = exp(log_sd_log_sao_dev);

  vector<Type> rec_pred = alpha * ssb / (beta + ssb);
  vector<Type> log_rec_pred = log(rec_pred);

  Type nll = 0.0;
  Type zero = 0.0;

  vector<Type> log_rec_resid = log_rec - log_rec_pred;
  vector<Type> log_rec_resid_std = log_rec_resid / sd_log_rec_me;

  // nll for log recruitment;
  nll -= dnorm(log_rec_resid, zero, sd_log_rec_me, true).sum();

  // nll for RW in lalpha_dev;
  nll -= dnorm(log_sao_dev(0), zero, sd_log_sao_dev, true);
  for (int i = 1; i < n; ++i)
  {
    nll -= dnorm(log_sao_dev(i), log_sao_dev(i - 1), sd_log_sao_dev, true);
  }

  REPORT(alpha);
  REPORT(beta);
  REPORT(sd_log_rec_me);

  REPORT(rec_pred);
  REPORT(log_rec_pred);
  REPORT(log_rec_resid);
  REPORT(log_rec_resid_std);

  vector<Type> log_sao_ts = log_sao + log_sao_dev;

  ADREPORT(alpha);
  ADREPORT(log_sao);
  ADREPORT(beta);
  ADREPORT(sd_log_rec_me);
  ADREPORT(sd_log_sao_dev);
  ADREPORT(log_sao_ts);

  return nll;
}
