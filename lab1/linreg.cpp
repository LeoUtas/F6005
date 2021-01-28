// Simple linear regression.
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x1);
  DATA_VECTOR(x2);

  PARAMETER(beta0);
  PARAMETER(beta1);
  PARAMETER(beta2);
  PARAMETER(logSigma);

  Type Sigma = exp(logSigma);
  Type Sigma2 = Sigma * Sigma;
  vector<Type> mu = beta0 + beta1 * x1 + beta2 * x2;

  Type nll = 0.0;

  nll = -sum(dnorm(Y, mu, Sigma, true));

  vector<Type> resid = Y - mu;
  vector<Type> std_resid = resid / Sigma;

  REPORT(Sigma);
  REPORT(mu);

  REPORT(resid);
  REPORT(std_resid);

  ADREPORT(mu);
  return nll;
}
