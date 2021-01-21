// Simple linear regression.
#include <TMB.hpp>
#include <math.h>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x1);
  DATA_VECTOR(x2);

  PARAMETER(beta_0);
  PARAMETER(beta_1);
  PARAMETER(beta_2);
  PARAMETER(logSigma);

  vector<Type> mu = beta_0 + beta_1 * x1 + beta_2 * x2;
  vector<Type> resid = Y - mu;

  Type nll = 0;
  Type Sigma = exp(logSigma);
  Type Sigma_2 = pow(2, Sigma);

  vector<Type> st_resid = resid / Sigma;

  nll = -dnorm(Y, mu, Sigma, true).sum();

  REPORT(resid);
  REPORT(st_resid);
  ADREPORT(Sigma_2);
  return nll;
}
