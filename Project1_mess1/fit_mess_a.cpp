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

  DATA_VECTOR(len_mid);
  DATA_VECTOR(len_border); // breaking border between length bins

  DATA_ARRAY(log_index);
  DATA_VECTOR(log_index2);
  DATA_VECTOR(log_index3);

  DATA_VECTOR(age);
  DATA_ARRAY(sf);

  DATA_MATRIX(weight);
  DATA_MATRIX(mat);

  DATA_IVECTOR(s1_iyear);
  DATA_IVECTOR(s1_ilen1);
  DATA_IVECTOR(s1_ilen2);

  DATA_IVECTOR(s2_iyear);
  DATA_IVECTOR(s2_ilen1);
  DATA_IVECTOR(s2_ilen2);

  DATA_IVECTOR(s3_iyear);
  DATA_IVECTOR(s3_ilen1);
  DATA_IVECTOR(s3_ilen2);

  DATA_IVECTOR(CL_iyear);
  DATA_IVECTOR(CL_ilen1);
  DATA_IVECTOR(CL_ilen2);
  DATA_VECTOR(log_catch);

  DATA_MATRIX(M);

  DATA_ARRAY(isurvey);

  array<Type> ni(isurvey);
  for (i = 0; i < isurvey; ++i)
  {
    ni(i) = index(i).size();
  }

  int nc = log_catch.size();
  Type zero = 0.0;
  Type one = 1.0;

  // *** declare parameters - fixed effects ***

  // for pla_func()
  PARAMETER(log_Linf);
  PARAMETER(log_vbk);
  PARAMETER(log_len_o);
  PARAMETER(log_cv_len);

  // for computing F
  PARAMETER(log_std_log_F);
  PARAMETER(logit_F_age);
  PARAMETER(logit_F_year);
  PARAMETER(log_F_main);

  // for cohort model
  PARAMETER(log_meanR);
  PARAMETER(log_std_log_R);
  PARAMETER(logit_log_R);
  PARAMETER_VECTOR(log_N0); // initial number at age (except for the first age), so size = A-1
  PARAMETER(log_std_pe);

  // for predicted survey index
  PARAMETER_VECTOR(log_std_index);
  PARAMETER_VECTOR(log_std_logq);

  // for predicted catch
  PARAMETER(log_std_catch);

  // *** declare parameters - random effects ***

  PARAMETER_VECTOR(log_Rec_dev);
  PARAMETER_ARRAY(log_F_dev);
  PARAMETER_ARRAY(pe);

  PARAMETER_ARRAY(logq);

  // *** transform parameters ***

  // for pla_func()
  Type Linf = exp(log_Linf);
  Type vbk = exp(log_vbk);
  Type cv_len = exp(log_cv_len);
  Type len_o = exp(log_len_o);
  Type po = len_o / Linf;

  // for computing F
  Type std_log_F = exp(log_std_log_F);
  Type phi_F_age = exp(logit_F_age) / (one + exp(logit_F_age));
  Type phi_F_year = exp(logit_F_year) / (one + exp(logit_F_year));

  // for cohort model
  Type std_log_R = exp(log_std_log_R);
  Type phi_logR = exp(logit_log_R) / (one + exp(logit_log_R));
  Type std_pe = exp(log_std_pe);

  // for predicted survey index
  Type std_index = exp(log_std_index);

  Type std_logq = exp(log_std_logq); // something can be revised here std_logq_vec ?

  // for predicted catch
  Type std_catch = exp(log_std_catch);

  using namespace density;

  Type nll = zero;
  int i, j;

  // *** compute length age key - proportion in each length bin, by age. mean_len and std_len not used in model;

  vector<Type> mean_len = Linf * (1 - (1 - po) * exp(-vbk * age));
  vector<Type> std_len = cv_len * mean_len;

  array<Type> age_s(isurvey);
  array<Type> pla_s(isurvey, L, A);

  for (i = 0; i < isurvey; ++i)
  {
    age_s(i) = age + sf(i);
    pla_s(i) = pla_func(L, A, Linf, po, vbk, cv_len, age_s(i), len_border);
  }

  // *** compute F ***

  matrix<Type> F(A, Y);
  matrix<Type> Z(A, Y);
  for (j = 0; j < Y; ++j)
  {
    F(0, j) = zero;
    F(1, j) = zero;
    for (i = 2; i < A; ++i)
    {
      F(i, j) = exp(log_F_main + log_F_dev(i - 2, j));
    }
  }
  Z = F + M;

  // log_F_dev
  array<Type> log_F_dev1(A - 3, Y); // A-3 because F_1 and F_2 = 0 , and F_A = F_A-1;
  for (i = 0; i < A - 3; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      log_F_dev1(i, j) = log_F_dev(i, j);
    }
  }

  // nll for F_dev ~ AR(1) for age and year
  nll += SCALE(SEPARABLE(AR1(phi_F_year), AR1(phi_F_age)), std_log_F)(log_F_dev1);
  // matrix exp() does not work

  // *** the cohort model ***

  vector<Type> log_Rec = log_Rec_dev + log_meanR;
  matrix<Type> log_N(A, Y);
  matrix<Type> N_matrix(A, Y);
  //initializing log numbers in  first year
  log_N(0, 0) = log_Rec(0); // log_Rec(0) ~ the number of the 1st age in the 1st year;
  for (i = 1; i < A; ++i)
  {
    log_N(i, 0) = log_N0(i - 1); // log_N) ~ the number of age i in the 1st year;
  }
  //compute log numbers at age
  for (j = 1; j < Y; ++j)
  {
    log_N(0, j) = log_Rec(j); // log_Rec(j) ~ the number of 1st age in year j;
    for (i = 1; i < A; ++i)
    {
      if (i < A - 1)
      {
        log_N(i, j) = log_N(i - 1, j - 1) - Z(i - 1, j - 1) + pe(i - 1, j - 1);
      }
      if (i == A - 1)
      {
        log_N(i, j) = log(exp(log_N(i - 1, j - 1) - Z(i - 1, j - 1)) + exp(log_N(i - 1, j) - Z(i - 1, j))) + pe(i - 1, j - 1);
      } //plus group
    }
  }
  N_matrix = exp(log_N.array());

  // nll for recruitment deviations in cohort model
  nll += SCALE(AR1(phi_logR), std_log_R)(log_Rec_dev);

  // nll for process error term
  for (int j = 0; j < Y - 1; ++j)
  {
    for (int i = 0; i < A - 1; ++i)
    {
      nll -= dnorm(pe(i, j), zero, std_pe, true);
    }
  }

  //compute numbers at time of survey (Ns);
  array<Type> Ns(isurvey, A, Y);

  for (i = 0; i < isurvey; ++i)
  {
    for (j = 0; j < A; ++j)
    {
      for (k = 0; k < Y; ++k)
      {
        Ns(i, j, k) = N_matrix(j, k) * exp(-sf(i) * Z(j, k));
      }
    }
  }

  // compute catch at age (CNA);
  matrix<Type> CNA(A, Y);

  for (i = 0; i < A; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      CNA(i, j) = N_matrix(i, j) * (one - exp(-Z(i, j))) * (F(i, j) / Z(i, j));
    }
  }

  //compute beginning of year pop num@len (NL), catch@len (CL), and survey num@len (NLs);
  matrix<Type> NL(L, Y);
  matrix<Type> pla_0 = pla_func(L, A, Linf, po, vbk, cv_len, age, len_border);
  NL = pla0 * N_matrix; // should use pla pla of age;

  matrix<Type> CL(L, Y);
  matrix<Type> pla_mid = pla_func(L, A, Linf, po, vbk, cv_len, age + .5, len_border);
  CL = pla_mid * CNA; // should use pla of (age + .5)

  array<Type> NLs(isurvey, L, Y);
  for (i = 0; i < isurvey; ++i)
  {
    NLs(i) = pla_s(i) * Ns(i);
  } // should use psf of age + sf;

  array<Type> log_NLs(isurvey);
  for (i = 0; i < isurvey; ++i)
  {
    log_NLs(i) = log(NLs(i).array());
  }

  //compute biomass, catch biomass and ssb @ len, and total biomass and ssb;
  matrix<Type> B_matrix = weight.array() * NL.array();      // stock biomass (in weight) at length for years;
  matrix<Type> CB_matrix = weight.array() * CL.array();     // catch (in weight) at length for years;
  matrix<Type> SSB_matrix = mat.array() * B_matrix.array(); // should revise SSB (use female only);
  vector<Type> biomass(Y);
  vector<Type> ssb(Y);
  biomass = B_matrix.colwise().sum();                     // stock biomass (sum of biomass at length) over years;
  ssb = SSB_matrix.colwise().sum();                       // sum of SSB at length over years;
  vector<Type> catch_biomass = CB_matrix.colwise().sum(); // sum of catch (in weight) at length over years;
  vector<Type> harvest_rate = catch_biomass / biomass;

  // ssb relative SSB to average, just FYI;
  vector<Type> rssb = Y * ssb / sum(ssb);

  // *** calculate predicted survey index and residuals ***
  array<Type> q(isurvey) = exp(logq(isurvey));
  array<Type> E_index(isurvey, ni(isurvey));

  for (i = 0; i < isruvey; ++i)
  {
    for (j = 0; j <= ni(isurvey); ++j)
    {
      E_index(i, j) = 0.0;
      for (k = s_ilen1(j); k <= s_ilen2(j); ++k)
      {
        E_index(i, j) += q(i, k) * NLs(k, s_iyear(j)); // should incorperate a pe in here;
      }
    }
  }

  array<Type> Elog_index(isurvey) = log(E_index(isurvey));
  array<Type> resid_index(isurvey) = log_index(isurvey) - Elog_index(isurvey);
  array<Type> std_resid_index(isurvey) = resid_index(isurvey) / std_index(isurvey);


  // nll -= dnorm(resid_index(isurvey), zero, std_index(isurvey), true).sum();

  // nll for catchability ~ logq 1
  for (i = 0; i < isurvey; ++i)
  {
    for (j = 3; j < 45; ++j)
    {
      nll -= dnorm(logq(i, j), logq(i, j - 1), std_logq(i), true);
    }
  }

  // *** calculate model predicted catch and residuals ***

  vector<Type> E_catch(nc);
  vector<Type> Elog_catch(nc);
  vector<Type> resid_catch(nc);
  vector<Type> std_resid_catch(nc);

  for (i = 0; i < nc; ++i)
  {
    E_catch(i) = 0.0;
    for (j = CL_ilen1(i); j <= CL_ilen2(i); ++j)
    {
      E_catch(i) += CL(j, CL_iyear(i)); // should incorperate a pe in here;
    }
  }
  Elog_catch = log(E_catch);
  resid_catch = log_catch - Elog_catch;
  std_resid_catch = resid_catch / std_catch;

  // nll for predicted catch
  nll -= dnorm(resid_catch, zero, std_catch, true).sum();

  REPORT(Linf);
  REPORT(vbk);
  REPORT(po);
  REPORT(len_o);
  REPORT(cv_len);

  REPORT(std_len);
  REPORT(mean_len);

  REPORT(pla_s);
  REPORT(pla_0);
  REPORT(pla_mid);

  REPORT(Z);
  REPORT(F);
  REPORT(std_log_F);

  REPORT(log_F_main);
  REPORT(log_F_dev);
  REPORT(phi_F_age);
  REPORT(phi_F_year);

  REPORT(log_N0);
  REPORT(N_matrix);

  REPORT(NL);

  REPORT(NLs);

  REPORT(Ns);

  REPORT(CNA);
  REPORT(CL);

  REPORT(log_meanR);
  REPORT(std_log_R);
  REPORT(log_Rec_dev);
  REPORT(log_Rec);
  REPORT(phi_logR);
  REPORT(pe);
  REPORT(std_pe);

  REPORT(B_matrix);
  REPORT(SSB_matrix);
  REPORT(CB_matrix);

  REPORT(biomass);
  REPORT(catch_biomass);
  REPORT(harvest_rate);
  REPORT(ssb);
  REPORT(rssb);

  REPORT(Elog_index);
  REPORT(E_index);
  REPORT(std_index);
  REPORT(resid_index);
  REPORT(std_resid_index);
  REPORT(logq);
  REPORT(std_logq);

  REPORT(E_catch);
  REPORT(Elog_catch);
  REPORT(resid_catch);
  REPORT(std_resid_catch);

  vector<Type> log_rssb = log(rssb); //log relative ssb;
  vector<Type> log_biomass = log(biomass);
  vector<Type> log_ssb = log(ssb);
  vector<Type> log_harvest_rate = log(harvest_rate);

  ADREPORT(log_biomass);
  ADREPORT(log_ssb);
  ADREPORT(log_rssb);
  ADREPORT(log_Rec);
  ADREPORT(log_harvest_rate);

  return nll;
}
