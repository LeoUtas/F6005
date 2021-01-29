#include <TMB.hpp>
#include <iostream>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(ssb);
  DATA_VECTOR(rec);
  DATA_VECTOR(lssb);
  DATA_VECTOR(lrec);

  PARAMETER(lalpha);
  PARAMETER(lbeta);
  PARAMETER(lsd_lrec_me);
  PARAMETER(lsd_lalpha_dev);

  PARAMETER_VECTOR(lalpha_dev);

  int n = ssb.size();

  Type beta = exp(lbeta);
  Type alpha_main = exp(lalpha);
  vector<Type> alpha = exp(lalpha + lalpha_dev);
  Type sd_lrec_me = exp(lsd_lrec_me);
  Type sd_lalpha_dev = exp(lsd_lalpha_dev);

  vector<Type> mu = alpha * ssb * exp(-beta * ssb);
  vector<Type> log_mu = log(mu);

  Type nll = 0.0;
  Type zero = 0.0;

  vector<Type> lrec_resid = lrec - log_mu;
  vector<Type> lrec_resid_std = lrec_resid / sd_lrec_me;

  // nll for log recruitment;
  nll -= sum(dnorm(lrec_resid, zero, sd_lrec_me, true));

  // nll for RW in lalpha_dev;
  nll -= dnorm(lalpha_dev(0), zero, sd_lalpha_dev, true);
  for (int i = 1; i < n; ++i)
  {
    nll -= dnorm(lalpha_dev(i), lalpha_dev(i - 1), sd_lalpha_dev, true);
  }

  vector<Type> lsao_ts = lalpha + lalpha_dev;

  REPORT(alpha_main);
  REPORT(beta);
  REPORT(alpha);
  REPORT(sd_lrec_me);

  REPORT(mu);
  REPORT(log_mu);
  REPORT(lrec_resid);
  REPORT(lrec_resid_std);

  ADREPORT(alpha_main);
  ADREPORT(beta);
  ADREPORT(sd_lrec_me);
  ADREPORT(sd_lalpha_dev);
  ADREPORT(lsao_ts);

  return nll;
}
