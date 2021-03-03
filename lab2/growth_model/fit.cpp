#include <TMB.hpp> 
#include <iostream>

template <class Type> 
Type pow(Type x,Type d){
  Type ret = exp(d*log(x));
  return ret;}

template<class Type>
Type objective_function<Type>::operator() ()
{   
  
  DATA_VECTOR(len);   
  DATA_VECTOR(age);   
  DATA_IVECTOR(igroup);
  DATA_IVECTOR(ispecies);
  DATA_IVECTOR(isex);    
  DATA_IVECTOR(istock);
  DATA_VECTOR(log_len);
  DATA_IVECTOR(p_igroup);  
  DATA_VECTOR(p_age);      
  DATA_IMATRIX(ui); 
  int nobs = len.size();
  int p_nobs = p_age.size();
  
//Fixed effects
   PARAMETER(log_Linfp); 
   PARAMETER(log_kp);    
   PARAMETER(log_deltap); 
   PARAMETER(log_Linfp_species); 
   PARAMETER(log_kp_species); 
   PARAMETER(log_Linfp_sex); 
   PARAMETER(log_kp_sex);  
    
   PARAMETER(log_std_log_Linf); 
   PARAMETER(log_std_log_k);   
   PARAMETER(log_std_log_delta);  
   PARAMETER(log_std_me_len);      
   PARAMETER(log_std_me_age);       
   PARAMETER_VECTOR(log_std_pe);     
   PARAMETER(logit_corr);  
   PARAMETER_VECTOR(log_alpha);     
   PARAMETER_VECTOR(log_beta);  
   
//Random Effects;  
   PARAMETER_VECTOR(log_Linf_stock); 
   PARAMETER_VECTOR(log_k_stock);   
   PARAMETER_VECTOR(log_delta_stock); 
   PARAMETER_VECTOR(pe);  
   PARAMETER_VECTOR(pe_pred); 
   PARAMETER_VECTOR(true_age);  
   int neffect = log_Linf_stock.size(); 
   
   Type zero = 0.0;
   Type one = 1.0;   
   Type two = 2.0; 
       
   Type std_log_Linf = exp(log_std_log_Linf);
   Type std_log_k = exp(log_std_log_k);           
   Type std_log_delta = exp(log_std_log_delta);
   Type deltap = exp(log_deltap);   
   Type std_me_len = exp(log_std_me_len);
   Type std_me_age = exp(log_std_me_age); 
   vector<Type> std_pe = exp(log_std_pe);
   Type corr = -exp(logit_corr)/(one + exp(logit_corr));    
   vector<Type> alpha = exp(log_alpha);  
   vector<Type> beta = exp(log_beta);    
   
   matrix<Type> covm(2,2);
   covm(0,0) =  exp(two*log_std_log_Linf);   
   covm(1,1) =  exp(two*log_std_log_k);
   covm(0,1) =  exp(log_std_log_k + log_std_log_Linf)*corr;
   covm(1,0) = covm(0,1);   
    
   vector<Type> pred(nobs);   
   vector<Type> log_pred(nobs);  
   vector<Type> log_pred_nope(nobs);
   vector<Type> resid(nobs);  
   vector<Type> std_resid(nobs);
   vector<Type> log_k(nobs); 
   vector<Type> log_Linf(nobs);  
   vector<Type> log_delta(nobs); 

//set up the sex and species differences;    
   vector<Type> log_Linfp_species_vec(2); 
   vector<Type> log_kp_species_vec(2);  
   vector<Type> log_Linfp_sex_vec(3); 
   vector<Type> log_kp_sex_vec(3);
   
   log_Linfp_species_vec(0) = zero; 
   log_Linfp_species_vec(1) = log_Linfp_species;  
   log_kp_species_vec(0) = zero; 
   log_kp_species_vec(1) = log_kp_species; 
   
   log_Linfp_sex_vec(0) = log_Linfp_sex/two;
   log_Linfp_sex_vec(1) = zero; 
   log_Linfp_sex_vec(2) = log_Linfp_sex; 
   log_kp_sex_vec(0) = log_kp_sex/two; 
   log_kp_sex_vec(1) = zero; 
   log_kp_sex_vec(2) = log_kp_sex;
           
   vector<Type> log_Linf_gp = log_Linfp + log_Linf_stock(ui.col(3)) + 
     log_Linfp_species_vec(ui.col(1)) + log_Linfp_sex_vec(ui.col(2));    
   vector<Type> log_k_gp = log_kp + log_k_stock(ui.col(3)) 
    + log_kp_species_vec(ui.col(1)) + log_kp_sex_vec(ui.col(2)); 
   vector<Type> log_delta_gp = log_deltap + log_delta_stock(ui.col(0));
          
   log_Linf = log_Linf_gp(igroup);
   log_k = log_k_gp(igroup);        
   log_delta = log_delta_gp(igroup);
     
   Type nll = zero;
    
   vector<Type> k = exp(log_k);   
   vector<Type> delta = exp(log_delta);
   log_pred = log_Linf + log(one - exp(-k*pow(true_age,delta))) + pe; 
   log_pred_nope = log_Linf + log(one - exp(-k*pow(true_age,delta))); 
   
   pred = exp(log_pred);
//   resid = len - pred;   
   resid = log(len/pred);
   std_resid = resid/std_me_len; 
   vector<Type> resid_age = log(age/true_age);
   vector<Type> std_resid_age = resid_age/std_me_age; 
   
// nll for len data;
   nll -= dnorm(resid, zero, std_me_len, true).sum();  
   
// nll for age data;
   nll -= dnorm(resid_age,zero, std_me_age, true).sum();
   nll -= dgamma(true_age,alpha(istock),beta(istock),true).sum();  
//   nll -= dnorm(log_true_age,log_alpha(istock),beta(istock),true).sum(); 

// nll for random delta stock effects;  
   nll -= dnorm(log_delta_stock, zero, std_log_delta, true).sum(); 

// nll for process errors;  
   nll -= dnorm(pe, zero, std_pe(istock), true).sum(); 
   vector<Type> sdg = std_pe(ui.col(3));
   nll -= dnorm(pe_pred, zero, sdg(p_igroup), true).sum();

// nll for random k,Linf effects; 
      
   vector<Type> vec_parm(2); 
      
   using namespace density;
   MVNORM_t<Type> neg_log_density(covm); 
   for(int i= 0; i < neffect; ++i){
     vec_parm(0) = log_Linf_stock(i); 
     vec_parm(1) = log_k_stock(i);     
     nll+=neg_log_density(vec_parm);
   }   
   
   vector<Type> p_log_k = log_k_gp(p_igroup);         
   vector<Type> p_log_Linf = log_Linf_gp(p_igroup); 
   vector<Type> p_log_delta = log_delta_gp(p_igroup);    
   vector<Type> p_k = exp(p_log_k);   
   vector<Type> p_delta = exp(p_log_delta);   
   Type da = Type(0.01);
   vector<Type> pd_age = p_age + da;
   vector<Type> p_log_pred = p_log_Linf + log(one - exp(-p_k*pow(pd_age,p_delta))) + pe_pred; 
   vector<Type> p_pred = exp(p_log_pred);
   
   REPORT(p_pred);
   REPORT(p_log_pred);
   REPORT(pred);
   REPORT(log_pred); 
   REPORT(log_pred_nope);
   REPORT(resid);
   REPORT(std_resid);
   REPORT(log_Linfp);   
   REPORT(log_kp);   
   REPORT(log_Linf_stock);   
   REPORT(log_k_stock);   
   REPORT(log_Linfp_species); 
   REPORT(log_kp_species); 
   REPORT(log_Linfp_sex); 
   REPORT(log_kp_sex);  
        
   REPORT(log_Linf);   
   REPORT(log_k);    
   REPORT(pe);  
   REPORT(true_age);   
   REPORT(resid_age);
   REPORT(std_resid_age);
      
   REPORT(std_log_Linf); 
   REPORT(std_log_k);    
   REPORT(std_pe);      
   REPORT(std_me_len);    
   REPORT(std_me_age);  
   REPORT(corr);        
   REPORT(covm);
   
   REPORT(log_Linf_gp);   
   REPORT(log_k_gp);      
   REPORT(log_delta_gp);      
   REPORT(log_Linfp);
   REPORT(log_kp); 
   REPORT(log_Linf_stock);   
   REPORT(log_k_stock);    
   
   ADREPORT(log_Linf_gp);   
   ADREPORT(log_k_gp);       
   ADREPORT(log_delta_gp);     
   ADREPORT(log_Linfp);
   ADREPORT(log_kp);   
   ADREPORT(log_Linf_stock);   
   ADREPORT(log_k_stock);      
   ADREPORT(p_pred);
   ADREPORT(p_log_pred); 

   return nll;
}

