#include <TMB.hpp>
#include <math.h>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_IVECTOR(year);
  DATA_VECTOR(C);
  DATA_VECTOR(index);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iq);
  DATA_VECTOR(log_C);
  DATA_VECTOR(log_index);

  // DATA_SCALAR(log_r_pred);
  // DATA_SCALAR(sd_log_r);
  DATA_SCALAR(log_Po_pred);
  DATA_SCALAR(sd_log_Po);
  DATA_SCALAR(sd_logC);
  // DATA_IVECTOR(istd); // ???

  PARAMETER(log_r);
  PARAMETER(log_K);
  PARAMETER_VECTOR(log_q);
  PARAMETER(log_Po);
  PARAMETER(log_Ho);
  PARAMETER(log_sd_rw);
  PARAMETER_VECTOR(log_sd_log_index);
  PARAMETER(log_sd_pe);
  PARAMETER(logit_ar_pe);
  PARAMETER(log_sd_logC); // ???

  PARAMETER_VECTOR(log_pe);
  PARAMETER_VECTOR(log_H_dev);

  int n = year.size();
  int ni = index.size();
  int i;
  Type one = 1.0;
  Type half = 0.5;
  Type zero = 0.0;

  // sd_logC = exp(log_sd_logC);

  Type r = exp(log_r);
  Type K = exp(log_K);
  Type sd_rw = exp(log_sd_rw);
  vector<Type> sd_log_index = exp(log_sd_log_index); //???
  Type sd_pe = exp(log_sd_pe);
  Type ar_pe = exp(logit_ar_pe) / (one + exp(logit_ar_pe)); // ???
  vector<Type> pe = exp(log_pe);

  vector<Type> log_P(n); //log population biomass divided by K at start of the year;
  // vector<Type> log_P_midy(n); // log P at middle of year;
  vector<Type> P(n);
  // vector<Type> P_midy(n);
  vector<Type> log_B(n);
  vector<Type> log_H(n);
  vector<Type> H(n);
  vector<Type> log_index_pred(ni);
  vector<Type> log_C_pred(n); // model log catch;

  Type nll = 0;

  // prior nll for log_r;
  //  nll -= dnorm(log_r,E_log_r,sd_log_r,true); // ???

  // prior nll for log_Po;
  nll -= dnorm(log_Po, log_Po_pred, sd_log_Po, true);

  // log of production model;

  P(0) = exp(log_Po);
  log_H(0) = log_Ho;
  H(0) = exp(log_H(0));
  for (i = 1; i < n; i++)
  {
    log_H(i) = log_H(i - 1) + log_H_dev(i - 1);
    H(i) = exp(log_H(i));
    P(i) = (P(i - 1) + r * P(i - 1) * (one - P(i - 1)) - H(i - 1) * P(i - 1)) * pe(i - 1);
  }

  log_P = log(P);
  log_B = log_K + log_P;
  log_C_pred = log_B + log_H;

  // for (i = 0; i < n - 1; i++)
  // {
  //   P_midy(i) = half * (P(i) + P(i + 1));
  // }

  // int ln = n - 1;
  // Type Pnp1 = (P(ln) + r * P(ln) * (one - P(ln)) - H(ln) * P(ln)) * pe(ln);
  // P_midy(ln) = half * (P(ln) + Pnp1);
  // log_P_midy = log(P_midy);

  // nll for index;

  log_index_pred = log_q(iq) + log_K + log_P(iyear);
  vector<Type> resid = log_index - log_index_pred;
  // vector<Type> sd_log_index_vec = sd_log_index(istd);
  vector<Type> std_resid = resid / sd_log_index;

  nll -= dnorm(resid, zero, sd_log_index, true).sum();

  // nll for catch;

  vector<Type> resid_C = log_C - log_C_pred;
  nll -= dnorm(resid_C, zero, sd_logC, true).sum();

  // nll for random walk deviation in log_H;

  nll -= dnorm(log_H_dev, zero, sd_rw, true).sum();

  //   nll for log_pe process errors;

  i = 0;
  nll -= dnorm(log_pe(i), zero, sd_pe / sqrt(one - pow(ar_pe,2)), true);
  for (int i = 1; i < n; ++i)
  {
    nll -= dnorm(log_pe(i) - ar_pe * log_pe(i - 1), zero, sd_pe, true);
  }

  // report output;

  Type Hmsy = half * r;
  Type Bmsy = half * K;
  Type MSY = Hmsy * Bmsy;
  vector<Type> log_rB = log_B - log(Bmsy);
  vector<Type> log_rH = log_H - log(Hmsy);
  vector<Type> B = exp(log_B);
  vector<Type> C_pred = exp(log_C_pred);
  vector<Type> index_pred = exp(log_index_pred);

  REPORT(log_r);
  REPORT(log_K);
  REPORT(log_q);
  REPORT(log_Po);
  REPORT(log_Ho);
  REPORT(Hmsy);
  REPORT(Bmsy);
  REPORT(MSY);
  REPORT(log_B);
  REPORT(log_H);
  REPORT(B);
  REPORT(H);
  REPORT(log_rB);
  REPORT(log_rH);
  REPORT(log_C_pred);
  REPORT(resid_C);
  REPORT(log_index_pred);
  REPORT(index_pred);
  REPORT(resid);
  REPORT(std_resid);
  REPORT(log_pe);
  REPORT(log_H_dev);

  ADREPORT(Hmsy);
  ADREPORT(Bmsy);
  ADREPORT(MSY);
  ADREPORT(ar_pe);
  ADREPORT(log_rB);
  ADREPORT(log_rH);
  ADREPORT(log_B);
  ADREPORT(log_H);

  return nll;
}
