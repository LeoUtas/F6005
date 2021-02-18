#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{

  // input data;

  DATA_VECTOR(log_index);
  DATA_VECTOR(log_q);
  DATA_VECTOR(sf);
  DATA_ARRAY(weight);
  DATA_ARRAY(mat);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_VECTOR(m);

  // fixed parameters;

  PARAMETER(log_mean_Rec);
  PARAMETER(log_sd_log_Rec);
  PARAMETER(log_sd_index);
  PARAMETER(log_sd_log_f);
  PARAMETER_VECTOR(log_sd_log_s);
  PARAMETER(log_sd_pe);
  PARAMETER(logit_log_Rec);
  PARAMETER_VECTOR(log_No);

  // random parameters;

  PARAMETER_VECTOR(log_Rec_dev);
  PARAMETER_VECTOR(log_sp);
  PARAMETER_VECTOR(log_f);
  PARAMETER_ARRAY(pe);

  // parameter transformations;

  int Y = mat.cols();
  int A = mat.rows();
  int n = log_index.size();

  Type zero = 0;
  Type one = 1;

  Type sd_log_Rec = exp(log_sd_log_Rec);
  Type sd_index = exp(log_sd_index);
  Type sd_log_f = exp(log_sd_log_f);
  vector<Type> sd_log_s = exp(log_sd_log_s);
  Type sd_pe = exp(log_sd_pe);
  Type phi_log_Rec = exp(logit_log_Rec) / (one + exp(logit_log_Rec));

  // set up some useful stuff;

  using namespace density;

  Type nll = zero;

  // set log_s = log_sp but with age 6 (i.e. tmb age 6) log_s set at 0;

  vector<Type> log_s(A);
  int i, j;
  log_s(0) = Type(-5.0);
  for (i = 1; i <= 5; ++i)
  {
    log_s(i) = log_sp(i - 1);
  }
  log_s(6) = zero;
  for (i = 7; i < A; ++i)
  {
    log_s(i) = log_sp(i - 2);
  }

  //compute F = f_y * s_a;

  array<Type> log_F(A, Y);
  array<Type> F(A, Y);
  array<Type> Z(A, Y);

  for (i = 0; i < A; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      log_F(i, j) = log_f(j) + log_s(i);
      F(i, j) = exp(log_F(i, j));
      Z(i, j) = F(i, j) + m(i);
    }
  }

  // The cohort model;
  //initializing first year

  array<Type> log_N(A, Y);

  for (i = 1; i < A; ++i)
  {
    log_N(i, 0) = log_No(i - 1);
  }

  //compute numbers at age

  vector<Type> log_Rec = log_Rec_dev + log_mean_Rec;

  log_N(0, 0) = log_Rec(0); // number at the 1st age in the 1st year;
  for (j = 1; j < Y; ++j)
  {
    log_N(0, j) = log_Rec(j); // number at the 1st age in the jth year;
    for (i = 1; i < A; ++i)
    {
      log_N(i, j) = log_N(i - 1, j - 1) - Z(i - 1, j - 1) + pe(i, j); // cohort model with pe
    }
  }

  // recruitment ~ AR(1);

  nll += SCALE(AR1(phi_log_Rec), sd_log_Rec)(log_Rec_dev);

  // process error ~ iid;

  for (int j = 0; j < Y; ++j)
  {
    for (int i = 0; i < A; ++i)
    {
      nll -= dnorm(pe(i, j), zero, sd_pe, true);
    }
  }

  // compute biomass;

  array<Type> N_matrix(A, Y);
  N_matrix = exp(log_N);
  array<Type> B_matrix = weight * N_matrix;
  array<Type> SSB_matrix = mat * B_matrix; // ???
  vector<Type> biomass(Y);
  vector<Type> ssb(Y);
  for (j = 0; j < Y; ++j)
  {
    biomass(j) = sum(B_matrix.col(j));
    ssb(j) = sum(SSB_matrix.col(j));
  }

  // compute ssb/average;

  vector<Type> rssb = Y * ssb / sum(ssb);

  // compute index;

  vector<Type> log_index_pred(n);
  vector<Type> resid_index(n);
  vector<Type> sd_resid_index(n);

  for (i = 0; i < n; ++i)
  {
    log_index_pred(i) = log_q(iage(i)) + log_N(iage(i), iyear(i)) - sf(i) * Z(iage(i), iyear(i));
  }
  resid_index = log_index - log_index_pred;
  sd_resid_index = resid_index / sd_index;

  // compute the joint negative log-likelihood of the indices;

  // index;

  nll -= dnorm(resid_index, zero, sd_index, true).sum();

  // s ~ random walk;

  Type sigma_s;

  for (i = 1; i < A - 1; ++i)
  {
    if (i < 6)
    {
      sigma_s = sd_log_s(0);
    }
    if (i >= 6)
    {
      sigma_s = sd_log_s(1);
    }
    nll -= dnorm(log_s(i), log_s(i - 1), sigma_s, true);
  }

  // f ~ random walk;

  for (i = 1; i < Y; ++i)
  {
    nll -= dnorm(log_f(i), log_f(i - 1), sd_log_f, true);
  };

  vector<Type> totN(Y);
  vector<Type> tni(Y);
  vector<Type> aveF_46(Y);
  for (i = 1; i < Y; ++i)
  {
    totN(i) = zero;
    aveF_46(i) = zero;
  }
  for (int i = 4; i <= 6; ++i)
  {
    tni = N_matrix.transpose().col(i);
    aveF_46 += vector<Type>(F.transpose().col(i)) * tni;
    totN += tni;
  }
  aveF_46 = aveF_46 / totN;
  vector<Type> log_aveF_46 = log(aveF_46);

  vector<Type> aveF_69(Y);
  for (i = 1; i < Y; ++i)
  {
    totN(i) = zero;
    aveF_69(i) = zero;
  }
  for (int i = 6; i <= 9; ++i)
  {
    tni = N_matrix.transpose().col(i);
    aveF_69 += vector<Type>(F.transpose().col(i)) * tni;
    totN += tni;
  }
  aveF_69 = aveF_69 / totN;
  vector<Type> log_aveF_69 = log(aveF_69);

  REPORT(N_matrix);
  REPORT(B_matrix);
  REPORT(SSB_matrix);
  REPORT(Z);
  REPORT(F);
  REPORT(log_s);
  REPORT(log_f);
  REPORT(biomass);
  REPORT(ssb);
  REPORT(rssb);
  REPORT(aveF_46);
  REPORT(aveF_69);
  REPORT(log_index_pred);
  REPORT(resid_index);
  REPORT(sd_resid_index);
  REPORT(log_Rec);
  REPORT(log_Rec_dev);
  REPORT(log_f);
  REPORT(pe);

  vector<Type> log_rssb = log(rssb);
  vector<Type> log_biomass = log(biomass);
  vector<Type> log_ssb = log(ssb);

  ADREPORT(log_biomass);
  ADREPORT(log_ssb);
  ADREPORT(log_rssb);
  ADREPORT(log_aveF_46);
  ADREPORT(log_aveF_69);
  ADREPORT(log_Rec);

  return nll;
}
