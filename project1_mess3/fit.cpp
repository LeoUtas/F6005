#include <TMB.hpp>
// first there are two functions, pnorm_discrete and pla_func, that should not be modified;

template <class Type>
Type pnorm_discrete(Type upper, Type lower)
{

  Type zero = 0.0;
  Type one = 1.0;

  Type mid = 0.5 * (upper + lower);
  Type mid2 = mid * mid;
  Type mid4 = mid * mid * mid * mid;
  Type ub = upper - mid;
  Type lb = lower - mid;
  Type d1 = upper - lower;
  Type d3 = ub * ub * ub - lb * lb * lb;
  Type d5 = ub * ub * ub * ub * ub - lb * lb * lb * lb * lb;
  Type ts = d1 + (mid2 - 1) * d3 / 6 + (mid4 - 6 * mid2 + 3) * d5 / 120; // taylor series expansion

  Type ret = dnorm(mid, zero, one, false) * ts;

  return ret;
}

template <class Type>
matrix<Type> pla_func(int L, int A, Type Linf, Type po, Type vbk, Type cv_len,
                      vector<Type> age, vector<Type> len_border)
{

  matrix<Type> pla(L, A);

  // compute length age key - proportion in each length bin, by age;
  for (int j = 0; j < A; ++j)
  {
    Type ml = Linf * (1 - (1 - po) * exp(-vbk * age(j)));
    Type sl = ml * cv_len;
    vector<Type> len_border_std = (len_border - ml) / sl; // standardized border among length bins

    int expand = 10; // the expansion of lower and upper bins
    vector<Type> lower_bin(expand);
    Type Ub = len_border_std(0);
    Type Lb = Ub - 1.0;
    lower_bin(0) = pnorm_discrete(Ub, Lb);
    for (int k = 1; k < expand; ++k)
    {
      Ub = Ub - 1.0;
      Lb = Ub - 1.0;
      lower_bin(k) = pnorm_discrete(Ub, Lb);
    }
    pla(0, j) = lower_bin.sum(); // first length bin

    for (int i = 1; i < (L - 1); ++i)
    {
      pla(i, j) = pnorm_discrete(len_border_std(i), len_border_std(i - 1));
    }

    vector<Type> upper_bin(expand);
    Lb = len_border_std(L - 2);
    Ub = Lb + 1.0;
    upper_bin(0) = pnorm_discrete(Ub, Lb);
    for (int k = 1; k < expand; ++k)
    {
      Lb = Lb + 1.0;
      Ub = Lb + 1.0;
      upper_bin(k) = pnorm_discrete(Ub, Lb);
    }
    pla(L - 1, j) = upper_bin.sum(); // last length bin
  }

  return pla;
}

// ########## the rest is the stock assessment model ##########################;
template <class Type>
Type objective_function<Type>::operator()()
{

  // input data;
  DATA_INTEGER(A);
  DATA_INTEGER(Y);
  DATA_INTEGER(L);
  DATA_VECTOR(age);
  DATA_VECTOR(len_mid);
  DATA_VECTOR(len_border); // breaking border between length bins

  DATA_VECTOR(log_index);
  DATA_VECTOR(sf);
  DATA_IVECTOR(is); //indicator for survey;
  DATA_MATRIX(weight);
  DATA_MATRIX(mat);
  DATA_VECTOR(F_ratio);
  DATA_VECTOR(juv_ratio);
  DATA_IVECTOR(SL_iyear);
  DATA_IVECTOR(SL_ilen1);
  DATA_IVECTOR(SL_ilen2);
  DATA_IVECTOR(CL_iyear);
  DATA_IVECTOR(CL_ilen1);
  DATA_IVECTOR(CL_ilen2);
  DATA_VECTOR(log_catch);
  DATA_MATRIX(M);

  int n = log_index.size();
  int nc = log_catch.size();
  int ns = sf.size(); //I assume a different sf for each survey so ns is number of survey indices;
  Type zero = 0.0;
  Type one = 1.0;
  Type half = 0.5;

  //define parameters;
  PARAMETER(log_meanR);
  PARAMETER(log_std_log_R);
  PARAMETER_VECTOR(log_std_index);
  PARAMETER(log_std_catch);
  PARAMETER_VECTOR(log_std_logq);
  PARAMETER(log_std_pe);
  PARAMETER(logit_log_R);
  PARAMETER(log_Linf);
  PARAMETER(log_vbk);
  PARAMETER(log_len_o);
  PARAMETER(log_cv_len);
  PARAMETER(log_std_log_F);
  PARAMETER(logit_F_age);
  PARAMETER(logit_F_year);
  PARAMETER_ARRAY(log_F_main);
  PARAMETER_VECTOR(log_N0); // initial number at age (except for the first age), so size = A-1

  // Random Effects
  PARAMETER_VECTOR(log_Rec_dev);
  PARAMETER_ARRAY(log_F_dev);
  PARAMETER_ARRAY(logq);
  PARAMETER_ARRAY(pe);

  // transform parameters;
  Type std_log_R = exp(log_std_log_R);
  vector<Type> std_index = exp(log_std_index);
  Type std_catch = exp(log_std_catch);
  vector<Type> std_logq = exp(log_std_logq);
  Type std_pe = exp(log_std_pe);
  Type std_log_F = exp(log_std_log_F);
  Type phi_logR = exp(logit_log_R) / (one + exp(logit_log_R));
  Type phi_F_age = exp(logit_F_age) / (one + exp(logit_F_age));
  Type phi_F_year = exp(logit_F_year) / (one + exp(logit_F_year));
  Type Linf = exp(log_Linf);
  Type vbk = exp(log_vbk);
  Type cv_len = exp(log_cv_len);
  Type len_o = exp(log_len_o);
  Type po = len_o / Linf;

  //make some containers;
  matrix<Type> log_N(A, Y);
  matrix<Type> Z(A, Y);
  matrix<Type> F(A, Y);
  matrix<Type> NL(L, Y);
  matrix<Type> CL(L, Y);
  array<Type> NLs(ns, L, Y);
  array<Type> log_NLs(ns, L, Y);
  matrix<Type> N_matrix(A, Y);
  array<Type> Ns(ns, A, Y);
  matrix<Type> CNA(A, Y);
  vector<Type> biomass(Y);
  vector<Type> mat_vec(Y);
  vector<Type> ssb(Y);
  vector<Type> juv(Y);
  vector<Type> Elog_index(n);
  vector<Type> E_index(n);
  vector<Type> resid_index(n);
  vector<Type> std_resid_index(n);
  vector<Type> E_catch(nc);
  vector<Type> Elog_catch(nc);
  vector<Type> resid_catch(nc);
  vector<Type> std_resid_catch(nc);

  using namespace density;

  Type nll = zero;
  int i, j, k;

  // *** Compute length age key - proportion ***;

  vector<Type> mean_len = Linf * (1 - (1 - po) * exp(-vbk * age));
  vector<Type> std_len = cv_len * mean_len;
  matrix<Type> pla = pla_func(L, A, Linf, po, vbk, cv_len, age, len_border); // pla for No;
  vector<Type> agep = age + half;
  matrix<Type> plac = pla_func(L, A, Linf, po, vbk, cv_len, agep, len_border); // pla for catch;
  array<Type> plas(ns, L, A);                                                  // pla for each sf;
  for (i = 0; i < ns; ++i)
  {
    vector<Type> agep = age + sf(i);
    matrix<Type> temp = pla_func(L, A, Linf, po, vbk, cv_len, agep, len_border);
    for (j = 0; j < L; ++j)
    {
      for (k = 0; k < A; ++k)
      {
        plas(i, j, k) = temp(j, k);
      }
    }
  }

  // *** Compute F ***;

  for (j = 0; j < Y; ++j)
  {
    F(0, j) = zero;
    F(1, j) = zero;
    for (i = 2; i < A; ++i)
    {
      F(i, j) = exp(log_F_main(i - 2, j) + log_F_dev(i - 2, j));
    }
  }
  Z = F + M;

  // log_F_dev nll
  array<Type> log_F_dev1(A - 3, Y); // A-3 because F_1 and F_2 = 0 , and F_A = F_A-1;
  for (i = 0; i < A - 3; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      log_F_dev1(i, j) = log_F_dev(i, j);
    }
  }

  nll += SCALE(SEPARABLE(AR1(phi_F_year), AR1(phi_F_age)), std_log_F)(log_F_dev1);

  // *** The cohort model ***;

  //A prior on the initial age distribution; steady state with Z=0.2, N_A = N_(A-1)exp(-Z);
  for (i = 1; i < A - 2; ++i)
  {
    nll -= dnorm(log_N0(i) - log_N0(i - 1), -Type(0.2), Type(0.5), true);
  }
  nll -= dnorm(log_N0(A - 2) - log_N0(A - 3), zero, Type(1.0), true);

  //initializing log numbers in  first year
  vector<Type> log_Rec = log_Rec_dev + log_meanR;
  log_N(0, 0) = log_Rec(0);
  for (i = 1; i < A; ++i)
  {
    log_N(i, 0) = log_N0(i - 1);
  }
  // compute log numbers at age
  for (j = 1; j < Y; ++j)
  {
    log_N(0, j) = log_Rec(j);
    for (i = 1; i < A; ++i)
    {
      if (i < A - 1)
      {
        log_N(i, j) = log_N(i - 1, j - 1) - Z(i - 1, j - 1) + pe(i - 1, j - 1);
      }
      if (i == A - 1)
      {
        log_N(i, j) = log(exp(log_N(i - 1, j - 1) - Z(i - 1, j - 1)) + exp(log_N(i, j - 1) - Z(i, j - 1))) + pe(i - 1, j - 1);
      } // plus group
    }
  }
  N_matrix = exp(log_N.array());

  // recruitment deviations
  nll += SCALE(AR1(phi_logR), std_log_R)(log_Rec_dev);

  // process error negative loglikelihood term
  for (int j = 0; j < Y - 1; ++j)
  {
    for (int i = 0; i < A - 1; ++i)
    {
      nll -= dnorm(pe(i, j), zero, std_pe, true);
    }
  }

  //compute numbers at time of catch at age (CNA)
  for (i = 0; i < A; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      CNA(i, j) = N_matrix(i, j) * (one - exp(-Z(i, j))) * (F(i, j) / Z(i, j));
    }
  }

  //compute numbers at time of survey (Ns)
  for (k = 0; k < ns; ++k)
  {
    for (i = 0; i < A; ++i)
    {
      for (j = 0; j < Y; ++j)
      {
        Ns(k, i, j) = N_matrix(i, j) * exp(-sf(k) * Z(i, j));
      }
    }
  }

  //compute beginning of year pop num@len (NL), catch@len (CL), and survey num@len (NLs);
  NL = pla * N_matrix;
  // for CL should use pla for mid year;
  CL = plac * CNA;
  // for Nls should use pla for time of year, sf;

  matrix<Type> temp1(L, A); //temp object to store pla for a survey
  matrix<Type> temp2(A, Y); //temp object to store Ns for a survey
  for (k = 0; k < ns; ++k)
  {
    //clunky code to extract matrices from the 3D arrays;
    for (i = 0; i < L; ++i)
    {
      for (j = 0; j < A; ++j)
      {
        temp1(i, j) = plas(k, i, j);
      }
    }
    for (i = 0; i < A; ++i)
    {
      for (j = 0; j < Y; ++j)
      {
        temp2(i, j) = Ns(k, i, j);
      }
    }

    matrix<Type> temp3 = temp1 * temp2;
    for (i = 0; i < L; ++i)
    {
      for (j = 0; j < Y; ++j)
      {
        NLs(k, i, j) = temp3(i, j);
        log_NLs(k, i, j) = log(NLs(k, i, j));
      }
    }
  }

  //compute biomass, catch biomass and mature biomass @ len, and total biomass and mature biomass;
  matrix<Type> B_matrix = weight.array() * NL.array();
  matrix<Type> CB_matrix = weight.array() * CL.array();
  matrix<Type> mat_matrix = mat.array() * B_matrix.array();
  biomass = B_matrix.colwise().sum();
  mat_vec = mat_matrix.colwise().sum();
  ssb = biomass * F_ratio;
  juv = biomass * juv_ratio;
  vector<Type> catch_biomass = CB_matrix.colwise().sum();
  vector<Type> harvest_rate = catch_biomass / biomass;

  // mature biomass and ssb relative to average over years, just FYI;
  vector<Type> rssb = Y * ssb / sum(ssb);
  vector<Type> rmat = Y * mat_vec / sum(mat_vec);

  // *** Compute predicted survey index & residuals ****;

  for (i = 0; i < n; ++i)
  {
    E_index(i) = 0.0;
    for (j = SL_ilen1(i); j <= SL_ilen2(i); ++j)
    {
      Type qi = exp(logq(is(i), j));
      E_index(i) += qi * NLs(is(i), j, SL_iyear(i));
    }
    Elog_index(i) = log(E_index(i));
    resid_index(i) = log_index(i) - Elog_index(i);
    std_resid_index(i) = resid_index(i) / std_index(is(i));
  }

  // nll of index
  nll -= dnorm(resid_index, zero, std_index(is), true).sum();

  // nll of logq
  for (k = 0; k < ns; ++k)
  {
    for (i = 3; i < 15; ++i)
    {
      nll -= dnorm(logq(k, i) - logq(k, i - 1), zero, Type(0.2), true);
    }
    for (i = 15; i < 30; ++i)
    {
      nll -= dnorm(logq(k, i), logq(k, i - 1), std_logq(k), true);
    }
    for (i = 30; i < 45; ++i)
    {
      nll -= dnorm(logq(k, i), logq(k, i - 1), Type(0.2), true);
    }
  }

  // // seem q and N are confounded; so I am adding medium priors on Fall RV q (len=20,22, 24) approx 1
  // nll -= dnorm(logq(0, 19), zero, Type(0.2), true);
  // nll -= dnorm(logq(0, 21), zero, Type(0.2), true);
  // nll -= dnorm(logq(0, 23), zero, Type(0.2), true);

  // *** Compute predicted catch and residuals *** ;
  for (i = 0; i < nc; ++i)
  {
    E_catch(i) = 0.0;
    for (j = CL_ilen1(i); j <= CL_ilen2(i); ++j)
    {
      E_catch(i) += CL(j, CL_iyear(i));
    }
  }
  Elog_catch = log(E_catch);
  resid_catch = log_catch - Elog_catch;
  std_resid_catch = resid_catch / std_catch;

  // nll of catch
  nll -= dnorm(resid_catch, zero, std_catch, true).sum();

  REPORT(pla);
  REPORT(plac);
  REPORT(plas);
  REPORT(mean_len);
  REPORT(log_len_o);
  REPORT(len_o);
  REPORT(cv_len);
  REPORT(Linf);
  REPORT(vbk);
  REPORT(po);
  REPORT(std_len);

  REPORT(log_meanR);
  REPORT(std_log_R);
  REPORT(logit_log_R);
  REPORT(phi_logR);
  REPORT(log_Rec_dev);
  REPORT(log_Rec);
  REPORT(pe);
  REPORT(std_pe);

  REPORT(log_N0);
  REPORT(N_matrix);
  REPORT(Ns);
  REPORT(NL);
  REPORT(NLs);
  REPORT(CNA);
  REPORT(CL);
  REPORT(B_matrix);
  REPORT(mat_matrix);
  REPORT(CB_matrix);
  REPORT(biomass);
  REPORT(catch_biomass);
  REPORT(ssb);
  REPORT(mat_vec);
  REPORT(rmat);

  REPORT(Z);
  REPORT(F);
  REPORT(std_log_F);
  REPORT(log_F_main);
  REPORT(phi_F_age);
  REPORT(phi_F_year);
  REPORT(logit_F_age);
  REPORT(logit_F_year);
  REPORT(log_F_dev);
  REPORT(harvest_rate);

  REPORT(Elog_index);
  REPORT(E_index);
  REPORT(resid_index);
  REPORT(std_resid_index);
  REPORT(std_index);
  REPORT(std_logq);
  REPORT(logq);

  REPORT(E_catch);
  REPORT(Elog_catch);
  REPORT(resid_catch);
  REPORT(std_resid_catch);

  vector<Type> log_rmat = log(rmat); //log relative mature biomass over years;
  vector<Type> log_biomass = log(biomass);
  vector<Type> log_mat_vec = log(mat_vec);
  vector<Type> log_ssb = log(ssb);
  vector<Type> log_juv = log(juv);
  vector<Type> log_harvest_rate = log(harvest_rate);

  ADREPORT(log_biomass);
  ADREPORT(log_mat_vec);
  ADREPORT(log_ssb);
  ADREPORT(log_juv);
  ADREPORT(log_rmat);
  ADREPORT(log_Rec);
  ADREPORT(log_harvest_rate);

  return nll;
}
