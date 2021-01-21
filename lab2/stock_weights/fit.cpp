#include <TMB.hpp> 
#include <iostream>

template<class Type>
  Type objective_function<Type>::operator() ()
{

  //input data; 
    DATA_VECTOR(x);  
    DATA_VECTOR(se);    
    DATA_IVECTOR(ia);
    DATA_IVECTOR(iy);
    DATA_IVECTOR(ic);   
    DATA_IVECTOR(iap);
    DATA_IVECTOR(iyp);
    DATA_IVECTOR(icp); 
  
    int n = x.size(); 
    int np = icp.size();  
    Type one = 1.0;
    Type zero = 0.0;  

  // parameter effects;  
    PARAMETER_VECTOR(age_eff);  
    PARAMETER_VECTOR(year_eff);
    PARAMETER_VECTOR(cohort_eff); 
    PARAMETER_VECTOR(log_std);  
    PARAMETER_VECTOR(logit_ar);  
    PARAMETER_ARRAY(dev);       
        
    vector<Type> ar = exp(logit_ar)/(one + exp(logit_ar));
    vector<Type> std = exp(log_std);
        
    Type ar_year_eff = ar(0); 
    Type ar_cohort_eff = ar(1);
    Type ar_dev_age = ar(2);  
    Type ar_dev_year = ar(3);
     
    Type std_year_eff = std(0);
    Type std_cohort_eff = std(1); 
    Type std_dev = std(2);
  
  //containers
     
    vector<Type> Ex(n); 
    vector<Type> resid(n); 
    vector<Type> std_resid(n);       
      
    //initialize the negative log likelihood
    Type nll = zero;   
    using namespace density;
     
    for(int i = 0;i < n;++i){
      Ex(i) = age_eff(ia(i)) + year_eff(iy(i)) + cohort_eff(ic(i)) + dev(iy(i),ia(i)); 
    }
    resid = x - Ex;
    std_resid = resid/se;
  
    
   //NEGATIVE LOGLIKELIHOODS
   //Index OBSERVATION MODEL
   nll -= dnorm(resid,zero,se,true).sum();
       
   //year effects
    nll += SCALE(AR1(ar_year_eff),std_year_eff)(year_eff);
       
   //cohort effects
    nll += SCALE(AR1(ar_cohort_eff),std_cohort_eff)(cohort_eff);
    
   // dev effects 
    nll += SCALE(SEPARABLE(AR1(ar_dev_age),AR1(ar_dev_year)),std_dev)(dev);
    
    
    vector<Type> log_pred_wt(np);
     
    for(int i = 0;i < np;++i){
      log_pred_wt(i) = age_eff(iap(i)) + year_eff(iyp(i)) + cohort_eff(icp(i)) + dev(iyp(i),iap(i)); 
    }
    
    
  	REPORT(age_eff); 
  	REPORT(year_eff); 
  	REPORT(cohort_eff); 
  	REPORT(dev);     
  	REPORT(Ex);    
  	REPORT(resid);     
  	REPORT(std_resid); 
    
  	REPORT(ar_year_eff);  
  	REPORT(ar_cohort_eff);  
  	REPORT(ar_dev_age);   
  	REPORT(ar_dev_year);  
   
  	REPORT(std_year_eff);  
  	REPORT(std_cohort_eff);  
  	REPORT(std_dev);   
                         
    ADREPORT(log_pred_wt); 
      
  return nll; 
  } 
