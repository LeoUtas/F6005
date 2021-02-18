#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{  
  DATA_VECTOR(ssb);  
  DATA_VECTOR(rec);  
  DATA_VECTOR(lssb);  
  DATA_VECTOR(lrec);
       
  PARAMETER(lsao);  
  PARAMETER(lbeta);  
  PARAMETER(lsd_lrec_me); 
  PARAMETER(lsd_lsao_dev);
  
  PARAMETER_VECTOR(lsao_dev);
  
  int n = ssb.size();   
   
  Type beta = exp(lbeta); 
  vector<Type> alpha = beta*exp(lsao+lsao_dev);    
  Type sd_lrec_me = exp(lsd_lrec_me); 
  Type sd_lsao_dev = exp(lsd_lsao_dev);  
      
  vector<Type> mu = alpha*ssb/(beta + ssb);
  vector<Type> log_mu = log(mu); 
     
  Type nll = 0.0; 
  Type zero=0.0;
  
  vector<Type> lrec_resid = lrec-log_mu;
  vector<Type> lrec_resid_std = lrec_resid/sd_lrec_me;  

// nll for log recruitment;  
  nll -= sum(dnorm(lrec_resid,zero,sd_lrec_me,true)); 
  
// nll for RW in lalpha_dev;   
  nll -= dnorm(lsao_dev(0),zero,sd_lsao_dev,true);
  for(int i = 1;i < n;++i){
    nll -= dnorm(lsao_dev(i),lsao_dev(i-1),sd_lsao_dev,true);
  }    
  
  REPORT(alpha);  
  REPORT(beta);  
  REPORT(sd_lrec_me); 
         
  REPORT(mu);             
  REPORT(log_mu);             
  REPORT(lrec_resid);          
  REPORT(lrec_resid_std);
  
  vector<Type> lsao_ts = lsao+lsao_dev;  
  
  ADREPORT(alpha);  
  ADREPORT(lsao);   
  ADREPORT(beta);  
  ADREPORT(sd_lrec_me); 
  ADREPORT(sd_lsao_dev); 
  ADREPORT(lsao_ts);     
 
  return nll;
}

