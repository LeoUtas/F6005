#include <TMB.hpp> 
#include <iostream>


template<class Type>
Type objective_function<Type>::operator() ()
{
 
  // input data;  
  DATA_VECTOR(log_index);
  DATA_VECTOR(logq);
  DATA_VECTOR(sf);      
  DATA_ARRAY(weight);
  DATA_ARRAY(mat);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage); 
  DATA_VECTOR(m);
  
  int Y = mat.cols(); 
  int A = mat.rows();  
  int n = log_index.size(); 
  Type zero = 0.0; 
  Type one = 1.0; 
  int i,j;

//define parameters;   
  PARAMETER(log_meanR);  
  PARAMETER(log_std_log_R);
  PARAMETER(log_std_index); 
  PARAMETER(log_meanF);  
  PARAMETER(log_std_log_F);
  PARAMETER_VECTOR(logit_log_F); 
  PARAMETER(log_std_pe);
  PARAMETER(logit_log_R); 
  PARAMETER_VECTOR(log_No);  
  PARAMETER_VECTOR(log_Rec_dev);
  PARAMETER_ARRAY(log_F_dev); 
  PARAMETER_ARRAY(pe);      
  
  Type std_log_R = exp(log_std_log_R); 
  Type std_index = exp(log_std_index); 
  Type std_log_F = exp(log_std_log_F); 
  Type std_pe = exp(log_std_pe); 
  Type phi_logR = exp(logit_log_R)/(one + exp(logit_log_R)); 
  vector<Type> phi_logF = exp(logit_log_F)/(one + exp(logit_log_F));
  Type phi_logF_year = phi_logF(0);           
  Type phi_logF_age = phi_logF(1);
  
  
  array<Type> log_N(A,Y); //log_N ages 0,...,A
  array<Type> log_F(A-1,Y); //log_F ages 1,...,A 
  array<Type> Z(A,Y);    
  array<Type> F(A,Y); 
  vector<Type> Elog_index(n); 
  vector<Type> resid_index(n);
  vector<Type> std_resid_index(n);  
   
  using namespace density;
  
  Type nll = zero;
  
  log_F = log_meanF + log_F_dev;
 
//compute F
  for(j = 0;j < Y;++j){
    F(0,j) = zero;    
    Z(0,j) = F(0,j) + m(0);
    for(i = 1;i < A;++i){
        F(i,j) = exp(log_F(i-1,j));
  	    Z(i,j) = F(i,j) + m(i);
     }
  }
 
  // The cohort model;
  //initializing first year  
  vector<Type> log_Rec = log_Rec_dev + log_meanR;
  log_N(0,0) = log_Rec(0);
  for(i = 1;i < A;++i){log_N(i,0) = log_No(i-1);} 
  //compute numbers at age
  for(j = 1;j < Y;++j){
    log_N(0,j) = log_Rec(j);  
    for(i = 1;i < A;++i){log_N(i,j) = log_N(i-1,j-1) - Z(i-1,j-1) + pe(i,j);}
  } 
  
  array<Type> N_matrix(A,Y);
  N_matrix = exp(log_N);
  array<Type> B_matrix = weight*N_matrix;  
  array<Type> SSB_matrix = mat*B_matrix; 
  vector<Type> biomass(Y);      
  vector<Type> ssb(Y);
  for(j = 0;j < Y;++j){
    biomass(j) = sum(B_matrix.col(j)); 
    ssb(j) = sum(SSB_matrix.col(j));
  }
  // ssb/average;
    vector<Type> rssb = Y*ssb/sum(ssb);

  //calculate index
  for(i = 0;i < n;++i){
  	Elog_index(i) = logq(iage(i)) + log_N(iage(i),iyear(i)) - sf(i)*Z(iage(i),iyear(i));
   }
 	 resid_index = log_index - Elog_index;
   std_resid_index = resid_index/std_index; 
 

  //negative log-likelihoods
  
  //index
   	nll -= dnorm(resid_index,zero,std_index,true).sum();
  
  //recruitment 
  nll += SCALE(AR1(phi_logR),std_log_R)(log_Rec_dev); 
  
  //nll for F devs
  nll += SCALE(SEPARABLE(AR1(phi_logF_year),AR1(phi_logF_age)),std_log_F)(log_F_dev);          

  //process error negative loglikelihood term
//  for(int j = 0;j < Y;++j){
//    for(int i = 0;i < A;++i){  
//      nll -= dnorm(pe(i,j),zero,std_pe,true);       
//    }
//  }  
    
  vector<Type> totN(Y);   
  vector<Type> tni(Y);   
  vector<Type> aveF_46(Y);
  for(i = 1;i < Y;++i){totN(i)=zero;aveF_46(i)=zero;}
  for(int i = 4;i <= 6;++i){
    tni = N_matrix.transpose().col(i);
    aveF_46 += vector<Type>(F.transpose().col(i))*tni;
    totN += tni;
  }
  aveF_46 = aveF_46/totN;
  vector<Type> log_aveF_46 = log(aveF_46); 
  
  vector<Type> aveF_69(Y);
  for(i = 1;i < Y;++i){totN(i)=zero;aveF_69(i)=zero;}
  for(int i = 6;i <= 9;++i){
    tni = N_matrix.transpose().col(i);
    aveF_69 += vector<Type>(F.transpose().col(i))*tni;
    totN += tni;
  }
  aveF_69 = aveF_69/totN;
  vector<Type> log_aveF_69 = log(aveF_69); 	
  
  REPORT(N_matrix); 
  REPORT(B_matrix);
  REPORT(SSB_matrix);
  REPORT(Z);          
  REPORT(F);
  REPORT(log_F_dev);  
  REPORT(biomass);
  REPORT(ssb);    
  REPORT(rssb);   
  REPORT(aveF_46); 
  REPORT(aveF_69);      
  REPORT(Elog_index); 
  REPORT(resid_index);  
  REPORT(std_resid_index); 
  REPORT(log_Rec);         
  REPORT(log_Rec_dev);
  REPORT(pe);         
  REPORT(phi_logF_year);
  REPORT(phi_logF_age);
  
  vector<Type> log_rssb = log(rssb);
  vector<Type> log_biomass = log(biomass);
  vector<Type> log_ssb = log(ssb);
  
  ADREPORT(log_biomass);
  ADREPORT(log_ssb);
  ADREPORT(log_rssb);
  ADREPORT(log_aveF_46); 
  ADREPORT(log_aveF_69);   
  ADREPORT(log_Rec);

  return nll;
}
