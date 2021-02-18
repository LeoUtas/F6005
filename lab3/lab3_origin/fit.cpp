#include <TMB.hpp> 

template<class Type>
Type objective_function<Type>::operator() ()
{  
  DATA_VECTOR(ssb);  
  DATA_VECTOR(rec);  
  DATA_VECTOR(lssb);  
  DATA_VECTOR(lrec);
       
  PARAMETER(lalpha);  
  PARAMETER(lbeta);  
  PARAMETER(lsd_lrec_me);
  PARAMETER(logit_ar_lrec_me);
  
  int n = ssb.size();   
  
  Type alpha = exp(lalpha);  
  Type beta = exp(lbeta);    
  Type sd_lrec_me = exp(lsd_lrec_me);  
  Type ar_lrec_me = exp(logit_ar_lrec_me)/(1.0 + exp(logit_ar_lrec_me));
  Type st_sd_lrec_me = sd_lrec_me/sqrt(1.0 - ar_lrec_me*ar_lrec_me);      
     
  Type nll = 0.0; 
  Type zero=0.0;
      
  vector<Type> mu = alpha*ssb/(beta + ssb);
  vector<Type> log_mu = log(mu); 
  
  vector<Type> lrec_resid = lrec-log_mu;
  vector<Type> lrec_resid_std = lrec_resid/st_sd_lrec_me;  

// nll for log recruitment;
  nll -= dnorm(lrec_resid(0),zero,st_sd_lrec_me,true);
// L34 does not work if 0.0 is used instead of zero??;   
  for(int i = 1;i < n;++i){
    nll -= dnorm(lrec_resid(i),ar_lrec_me*lrec_resid(i-1),sd_lrec_me,true);
  }
  
//  using namespace density;
//  nll += SCALE(AR1(ar_lrec_me),st_sd_lrec_me)(lrec_resid);    
  
  REPORT(alpha);  
  REPORT(beta);  
  REPORT(sd_lrec_me); 
  REPORT(ar_lrec_me); 
  REPORT(st_sd_lrec_me);
         
  REPORT(mu);             
  REPORT(log_mu);             
  REPORT(lrec_resid);          
  REPORT(lrec_resid_std); 
  
  ADREPORT(alpha);  
  ADREPORT(beta);  
  ADREPORT(sd_lrec_me); 
  ADREPORT(ar_lrec_me); 
  ADREPORT(st_sd_lrec_me);     
 
  return nll;
}

