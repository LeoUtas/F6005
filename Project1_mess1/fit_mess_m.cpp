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
  //DATA_VECTOR(len_mid);
  DATA_VECTOR(len_border); // breaking border between length bins

  DATA_VECTOR(log_index1);
  DATA_VECTOR(log_index2);
  DATA_VECTOR(log_index3);
  DATA_VECTOR(log_index4);

  DATA_SCALAR(sf1);
  DATA_SCALAR(sf2);
  DATA_SCALAR(sf3);
  DATA_SCALAR(sf4);

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

  DATA_IVECTOR(s4_iyear);
  DATA_IVECTOR(s4_ilen1);
  DATA_IVECTOR(s4_ilen2);

  DATA_IVECTOR(CL_iyear);
  DATA_IVECTOR(CL_ilen1);
  DATA_IVECTOR(CL_ilen2);
  DATA_VECTOR(log_catch);

  DATA_MATRIX(M);

  int ni1 = log_index1.size();
  int ni2 = log_index2.size();
  int ni3 = log_index3.size();
  int ni4 = log_index4.size();

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
  PARAMETER(log_std_index1);
  PARAMETER(log_std_index2);
  PARAMETER(log_std_index3);
  PARAMETER(log_std_index4);

  PARAMETER(log_std_logq1);
  PARAMETER(log_std_logq2);
  PARAMETER(log_std_logq3);
  PARAMETER(log_std_logq4);

  // for predicted catch
  PARAMETER(log_std_catch);

  // *** declare parameters - random effects ***

  PARAMETER_VECTOR(log_Rec_dev);
  PARAMETER_ARRAY(log_F_dev);
  PARAMETER_ARRAY(pe);

  PARAMETER_VECTOR(logq1);
  PARAMETER_VECTOR(logq2);
  PARAMETER_VECTOR(logq3);
  PARAMETER_VECTOR(logq4);

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
  Type std_index1 = exp(log_std_index1);
  Type std_index2 = exp(log_std_index2);
  Type std_index3 = exp(log_std_index3);
  Type std_index4 = exp(log_std_index4);

  Type std_logq1 = exp(log_std_logq1); // something can be revised here std_logq_vec ?
  Type std_logq2 = exp(log_std_logq2);
  Type std_logq3 = exp(log_std_logq3);
  Type std_logq4 = exp(log_std_logq4);

  // for predicted catch
  Type std_catch = exp(log_std_catch);

  // *** declare containers to store outputs ***

  // for computing F
  matrix<Type> F(A, Y);
  matrix<Type> Z(A, Y);

  // for cohort model
  matrix<Type> log_N(A, Y);
  matrix<Type> N_matrix(A, Y);

  matrix<Type> NL(L, Y);
  matrix<Type> Ns1(A, Y);
  matrix<Type> Ns2(A, Y);
  matrix<Type> Ns3(A, Y);
  matrix<Type> Ns4(A, Y);

  matrix<Type> CNA(A, Y);
  matrix<Type> CL(L, Y);

  matrix<Type> NLs1(L, Y);
  matrix<Type> NLs2(L, Y);
  matrix<Type> NLs3(L, Y);
  matrix<Type> NLs4(L, Y);

  vector<Type> biomass(Y);
  vector<Type> ssb(Y);

  // for predicted survey index
  vector<Type> Elog_index1(ni1);
  vector<Type> Elog_index2(ni2);
  vector<Type> Elog_index3(ni3);
  vector<Type> Elog_index4(ni4);

  vector<Type> E_index1(ni1);
  vector<Type> E_index2(ni2);
  vector<Type> E_index3(ni3);
  vector<Type> E_index4(ni4);

  vector<Type> resid_index1(ni1);
  vector<Type> resid_index2(ni2);
  vector<Type> resid_index3(ni3);
  vector<Type> resid_index4(ni4);

  vector<Type> std_resid_index1(ni1);
  vector<Type> std_resid_index2(ni2);
  vector<Type> std_resid_index3(ni3);
  vector<Type> std_resid_index4(ni4);

  // for predicted catch
  vector<Type> E_catch(nc);
  vector<Type> Elog_catch(nc);
  vector<Type> resid_catch(nc);
  vector<Type> std_resid_catch(nc);

  using namespace density;

  Type nll = zero;
  int i, j;

  // compute length age key - proportion in each length bin, by age. mean_len and std_len not used in model;

  vector<Type> mean_len = Linf * (1 - (1 - po) * exp(-vbk * age));
  vector<Type> std_len = cv_len * mean_len;

  vector<Type> age_s1 = age + sf1;
  vector<Type> age_s2 = age + sf2;
  vector<Type> age_s3 = age + sf3;
  vector<Type> age_s4 = age + sf4;

  matrix<Type> pla1 = pla_func(L, A, Linf, po, vbk, cv_len, age_s1, len_border);
  matrix<Type> pla2 = pla_func(L, A, Linf, po, vbk, cv_len, age_s2, len_border);
  matrix<Type> pla3 = pla_func(L, A, Linf, po, vbk, cv_len, age_s3, len_border);
  matrix<Type> pla4 = pla_func(L, A, Linf, po, vbk, cv_len, age_s4, len_border);

  // *** compute F ***

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

  //compute numbers at time of survey (Ns) and catch at age (CNA)
  for (i = 0; i < A; ++i)
  {
    for (j = 0; j < Y; ++j)
    {
      Ns1(i, j) = N_matrix(i, j) * exp(-sf1 * Z(i, j));
      Ns2(i, j) = N_matrix(i, j) * exp(-sf2 * Z(i, j));
      Ns3(i, j) = N_matrix(i, j) * exp(-sf3 * Z(i, j));
      Ns4(i, j) = N_matrix(i, j) * exp(-sf4 * Z(i, j));
      CNA(i, j) = N_matrix(i, j) * (one - exp(-Z(i, j))) * (F(i, j) / Z(i, j));
    }
  }

  //compute beginning of year pop num@len (NL), catch@len (CL), and survey num@len (NLs);
  matrix<Type> pla0 = pla_func(L, A, Linf, po, vbk, cv_len, age, len_border);
  NL = pla0 * N_matrix;

  vector<Type> age_mid = age + .5;
  matrix<Type> pla_mid = pla_func(L, A, Linf, po, vbk, cv_len, age_mid, len_border);
  CL = pla_mid * CNA; // use pla of (age + .5)
  // for Nls should use pla for time of year, sf;
  NLs1 = pla1 * Ns1;
  NLs2 = pla2 * Ns2;
  NLs3 = pla3 * Ns3;
  NLs4 = pla4 * Ns4; // use pla of (age + sf)

  matrix<Type> log_NLs1 = log(NLs1.array());
  matrix<Type> log_NLs2 = log(NLs2.array());
  matrix<Type> log_NLs3 = log(NLs3.array());
  matrix<Type> log_NLs4 = log(NLs4.array());

  //compute biomass, catch biomass and ssb @ len, and total biomass and ssb;
  matrix<Type> B_matrix = weight.array() * NL.array();
  matrix<Type> CB_matrix = weight.array() * CL.array();
  matrix<Type> SSB_matrix = mat.array() * B_matrix.array(); // should revise SSB (use female only);
  biomass = B_matrix.colwise().sum();
  ssb = SSB_matrix.colwise().sum(); // should revise SSB (use female only);
  vector<Type> catch_biomass = CB_matrix.colwise().sum();
  vector<Type> harvest_rate = catch_biomass / biomass;

  // ssb relative SSB to average, just FYI;
  vector<Type> rssb = Y * ssb / sum(ssb);

  // *** calculate predicted survey index and residuals ***

  vector<Type> q1 = exp(logq1);
  vector<Type> q2 = exp(logq2);
  vector<Type> q3 = exp(logq3);
  vector<Type> q4 = exp(logq4);

  // survey index 1
  for (i = 0; i < ni1; ++i)
  {
    E_index1(i) = 0;
    for (j = s1_ilen1(i); j <= s1_ilen2(i); ++j)
    {
      E_index1(i) += q1(j) * NLs1(j, s1_iyear(i)); // a possible pe in here;
    }
  }

  // survey index 2
  for (i = 0; i < ni2; ++i)
  {
    E_index2(i) = 0;
    for (j = s2_ilen1(i); j <= s2_ilen2(i); ++j)
    {
      E_index2(i) += q2(j) * NLs2(j, s2_iyear(i)); // a possible pe in here;
    }
  }

  // survey index 3
  for (i = 0; i < ni3; ++i)
  {
    E_index3(i) = 0;
    for (j = s3_ilen1(i); j <= s3_ilen2(i); ++j)
    {
      E_index3(i) += q3(j) * NLs3(j, s3_iyear(i)); // a possible pe in here;;
    }
  }

  // survey index 4
  for (i = 0; i < ni4; ++i)
  {
    E_index4(i) = 0;
    for (j = s4_ilen1(i); j <= s4_ilen2(i); ++j)
    {
      E_index4(i) += q4(j) * NLs4(j, s4_iyear(i)); // a possible pe in here;;
    }
  }

  // survey index 1
  Elog_index1 = log(E_index1);
  resid_index1 = log_index1 - Elog_index1;
  std_resid_index1 = resid_index1 / std_index1;

  // nll  for survey index 1
  nll -= dnorm(resid_index1, zero, std_index1, true).sum();

  // survey index 2
  Elog_index2 = log(E_index2);
  resid_index2 = log_index2 - Elog_index2;
  std_resid_index2 = resid_index2 / std_index2;

  // nll  for survey index 2
  nll -= dnorm(resid_index2, zero, std_index2, true).sum();

  // survey index 3
  Elog_index3 = log(E_index3);
  resid_index3 = log_index3 - Elog_index3;
  std_resid_index3 = resid_index3 / std_index3;

  // nll  for survey index 3
  nll -= dnorm(resid_index3, zero, std_index3, true).sum();

  // survey index 4
  Elog_index4 = log(E_index4);
  resid_index4 = log_index4 - Elog_index4;
  std_resid_index4 = resid_index4 / std_index4;

  // nll  for survey index 4
  nll -= dnorm(resid_index4, zero, std_index4, true).sum();

  // nll for catchability ~ logq 1
  for (i = 3; i < 26; ++i)
  {
    nll -= dnorm(logq1(i), logq1(i - 1), std_logq1, true);
  }

  // nll for catchability ~ logq 2
  for (i = 3; i < 26; ++i)
  {
    nll -= dnorm(logq2(i), logq2(i - 1), std_logq2, true);
  }

  // nll for catchability ~ logq 3
  for (i = 5; i < 26; ++i)
  {
    nll -= dnorm(logq3(i), logq3(i - 1), std_logq3, true);
  }

  // nll for catchability ~ logq 4
  for (i = 5; i < 26; ++i)
  {
    nll -= dnorm(logq4(i), logq4(i - 1), std_logq4, true);
  }

  // *** calculate model predicted catch and residuals ***

  for (i = 0; i < nc; ++i)
  {
    E_catch(i) = 0;
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

  REPORT(pla1);
  REPORT(pla2);
  REPORT(pla3);
  REPORT(pla4);

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

  REPORT(NLs1);
  REPORT(NLs2);
  REPORT(NLs3);
  REPORT(NLs4);

  REPORT(Ns1);
  REPORT(Ns2);
  REPORT(Ns3);
  REPORT(Ns4);

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

  REPORT(Elog_index1);
  REPORT(Elog_index2);
  REPORT(Elog_index3);
  REPORT(Elog_index4);

  REPORT(E_index1);
  REPORT(E_index2);
  REPORT(E_index3);
  REPORT(E_index4);

  REPORT(std_index1);
  REPORT(std_index2);
  REPORT(std_index3);
  REPORT(std_index4);

  REPORT(resid_index1);
  REPORT(resid_index2);
  REPORT(resid_index3);
  REPORT(resid_index4);

  REPORT(std_resid_index1);
  REPORT(std_resid_index2);
  REPORT(std_resid_index3);
  REPORT(std_resid_index4);

  REPORT(logq1);
  REPORT(logq2);
  REPORT(logq3);
  REPORT(logq4);

  REPORT(std_logq1);
  REPORT(std_logq2);
  REPORT(std_logq3);
  REPORT(std_logq4);

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
