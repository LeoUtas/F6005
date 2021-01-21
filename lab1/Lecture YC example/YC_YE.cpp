#include <TMB.hpp>
#include <iostream>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(iAE);
  DATA_IVECTOR(iCE);
  DATA_IVECTOR(iYE)

  // fixed effects;
  PARAMETER_VECTOR(par_AE);
  PARAMETER_VECTOR(par_CE);
  PARAMETER(log_std_index);
  PARAMETER(log_std_YE);

  // random effects;
  PARAMETER_VECTOR(par_YE);

  //  int n = y.size();
  //  int npar_ye = par_YE.size();
  int npar_ae = par_AE.size();
  int i;
  //  int i,j;

  Type std_index = exp(log_std_index);
  Type std_YE = exp(log_std_YE);
  vector<Type> par_AE_all(npar_ae + 1);

  for (i = 0; i < npar_ae; i++)
  {
    par_AE_all(i + 1) = par_AE(i);
  }
  par_AE_all(0) = 0.0;

  Type nll = 0;
  vector<Type> Ey = par_CE(iCE) + par_AE_all(iAE) + par_YE(iYE);

  // nll for observed index;
  nll -= sum(dnorm(y, Ey, std_index, true));

  //   nll for YE effect;
  nll -= sum(dnorm(par_YE, 0.0, std_YE, true));

  vector<Type> resid = y - Ey;
  vector<Type> std_resid = resid / std_index;

  REPORT(Ey);
  REPORT(resid);
  REPORT(std_resid);
  REPORT(par_AE);
  REPORT(par_AE_all);
  REPORT(par_CE);
  REPORT(par_YE);
  REPORT(nll);

  vector<Type> CE = exp(par_CE);
  Type mean_CE = sum(CE) / CE.size();
  vector<Type> CE_dev = log(CE / mean_CE);

  ADREPORT(Ey);
  ADREPORT(CE_dev);

  return nll;
}
