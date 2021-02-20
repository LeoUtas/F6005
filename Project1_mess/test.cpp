#include <TMB.hpp> 
// first there are two functions, pnorm_discrete and pla_func, that should not be modified;

template<class Type>
Type pnorm_discrete(Type upper, Type lower){
	
	Type zero = 0.0;
	Type one = 1.0;
	
	Type mid = 0.5 * (upper + lower);
	Type mid2 = mid*mid;
	Type mid4 = mid*mid*mid*mid;
	Type ub = upper - mid;
	Type lb = lower - mid;
	Type d1 = upper - lower;
	Type d3 = ub*ub*ub - lb*lb*lb;
	Type d5 = ub*ub*ub*ub*ub - lb*lb*lb*lb*lb;
	Type ts = d1 + (mid2-1)*d3/6 + (mid4-6*mid2+3)*d5/120; // taylor series expansion
	
	Type ret = dnorm(mid,zero,one,false) * ts;
	
	return ret;
}


template<class Type>
matrix<Type> pla_func(int L, int A, Type Linf, Type po, Type vbk, Type cv_len,
vector<Type> age, vector<Type> len_border){
  
  matrix<Type> pla(L,A);
  
  // compute length age key - proportion in each length bin, by age;
	for(int j = 0;j < A;++j){
		Type ml = Linf*(1-(1-po)*exp(-vbk*age(j)));
		Type sl = ml*cv_len;
		vector<Type> len_border_std = (len_border-ml)/sl; // standardized border among length bins
		
		int expand = 10; // the expansion of lower and upper bins
		vector<Type> lower_bin(expand);
		Type Ub = len_border_std(0);
		Type Lb = Ub - 1.0;
		lower_bin(0) = pnorm_discrete(Ub,Lb);
		for(int k=1; k<expand; ++k){
			Ub = Ub - 1.0;
			Lb = Ub - 1.0;
			lower_bin(k) = pnorm_discrete(Ub,Lb);
		}
		pla(0,j) = lower_bin.sum(); // first length bin
		
		for(int i=1; i< (L-1);++i){
			pla(i,j) = pnorm_discrete(len_border_std(i),len_border_std(i-1));
		}
		
		vector<Type> upper_bin(expand);
		Lb = len_border_std(L-2);
		Ub = Lb + 1.0;
		upper_bin(0) = pnorm_discrete(Ub,Lb);
		for(int k=1; k<expand; ++k){
			Lb = Lb + 1.0;
			Ub = Lb + 1.0;
			upper_bin(k) = pnorm_discrete(Ub,Lb);
		}
		pla(L-1,j) = upper_bin.sum(); // last length bin		
	}
	
	return pla;
}

// ########## the rest is the stock assessment model ##########################;
template<class Type>
Type objective_function<Type>::operator() ()
{
 
  // input data;
  DATA_INTEGER(A); 
  DATA_INTEGER(Y); 
  DATA_INTEGER(L); 
  DATA_VECTOR(age); 
  DATA_VECTOR(len_mid); 
  DATA_VECTOR(len_border); // breaking border between length bins  
  
  DATA_VECTOR(log_index);
  DATA_SCALAR(sf);      
  DATA_MATRIX(weight);
  DATA_MATRIX(mat);
  DATA_IVECTOR(FRV_iyear);
  DATA_IVECTOR(FRV_ilen1);   
  DATA_IVECTOR(FRV_ilen2);  
  DATA_IVECTOR(CL_iyear);
  DATA_IVECTOR(CL_ilen1);   
  DATA_IVECTOR(CL_ilen2); 
  DATA_VECTOR(log_catch);
  DATA_MATRIX(M);  
   
  int n = log_index.size(); 
  int nc = log_catch.size(); 
  Type zero = 0.0; 
  Type one = 1.0; 

//define parameters;   
  PARAMETER(log_meanR);  
  PARAMETER(log_std_log_R);
  PARAMETER(log_std_index);  
  PARAMETER(log_std_catch); 
  PARAMETER(log_std_logq);   
  PARAMETER(log_std_pe);
  PARAMETER(logit_log_R); 
  PARAMETER(log_Linf);
  PARAMETER(log_vbk);  
  PARAMETER(log_len_o);
  PARAMETER(log_cv_len); 
  PARAMETER(log_std_log_F); 
  PARAMETER(logit_F_age); 
  PARAMETER(logit_F_year);  
  PARAMETER(log_F_main);  
  PARAMETER_VECTOR(log_N0); // initial number at age (except for the first age), so size = A-1 
// Random Effects  
  PARAMETER_VECTOR(log_Rec_dev);   
  PARAMETER_ARRAY(log_F_dev);  
  PARAMETER_VECTOR(logq); 
  PARAMETER_ARRAY(pe);      

// transform parameters;  
  Type std_log_R = exp(log_std_log_R); 
  Type std_index = exp(log_std_index);  
  Type std_catch = exp(log_std_catch);   
  Type std_logq = exp(log_std_logq); 
  Type std_pe = exp(log_std_pe); 
  Type std_log_F = exp(log_std_log_F); 
  Type phi_logR = exp(logit_log_R)/(one + exp(logit_log_R)); 
  Type phi_F_age = exp(logit_F_age)/(one + exp(logit_F_age));
  Type phi_F_year = exp(logit_F_year)/(one + exp(logit_F_year));  
  Type Linf = exp(log_Linf);
  Type vbk = exp(log_vbk);  
  Type cv_len = exp(log_cv_len); 
  Type len_o = exp(log_len_o);
  Type po = len_o/Linf;

//make some containers;  
  matrix<Type> log_N(A,Y); 
  matrix<Type> Z(A,Y);    
  matrix<Type> F(A,Y); 
  matrix<Type> NL(L,Y); 
  matrix<Type> CL(L,Y); 
  matrix<Type> NLs(L,Y); 
  matrix<Type> N_matrix(A,Y);
  matrix<Type> Ns(A,Y);  
  matrix<Type> CNA(A,Y);  
  vector<Type> biomass(Y);      
  vector<Type> ssb(Y);
  vector<Type> Elog_index(n); 
  vector<Type> E_index(n); 
  vector<Type> resid_index(n);
  vector<Type> std_resid_index(n); 
  vector<Type> E_catch(nc);      
  vector<Type> Elog_catch(nc); 
  vector<Type> resid_catch(nc);
  vector<Type> std_resid_catch(nc);
  
  using namespace density; 
  
  Type nll = zero;
  int i,j;
  
  // compute length age key - proportion in each length bin, by age. mean_len and std_len not used in model;
  vector<Type> mean_len = Linf*(1-(1-po)*exp(-vbk*age));
  vector<Type> std_len = cv_len*mean_len;  
  matrix<Type> pla = pla_func(L,A,Linf,po,vbk,cv_len,age,len_border);
 
//compute F; matrix exp() does not work 
  for(j = 0;j < Y;++j){ 
    F(0,j) = zero;       
    F(1,j) = zero;
    for(i = 2;i < A;++i){
        F(i,j) = exp(log_F_main + log_F_dev(i-2,j));
  }}
  Z = F + M;
  
  vector<Type> log_Rec = log_Rec_dev + log_meanR;
 
  // The cohort model;
  //initializing log numbers in  first year 
  log_N(0,0) = log_Rec(0);
  for(i = 1;i < A;++i){log_N(i,0) = log_N0(i-1);}
  //compute log numbers at age
  for(j = 1;j < Y;++j){
    log_N(0,j) = log_Rec(j);  
    for(i = 1;i < A;++i){
      if(i<A-1){log_N(i,j) = log_N(i-1,j-1) - Z(i-1,j-1) + pe(i-1,j-1);}     
      if(i==A-1){log_N(i,j) = log(exp(log_N(i-1,j-1)-Z(i-1,j-1))+exp(log_N(i-1,j)-Z(i-1,j)))+pe(i-1,j-1);} //plus group
    } 
  }   
  N_matrix = exp(log_N.array());
  
//compute numbers at time of survey (Ns) and catch at age (CNA)
  for(i = 0;i < A;++i){  
    for(j = 0;j < Y;++j){
      Ns(i,j) = N_matrix(i,j)*exp(-sf*Z(i,j)); 
	  	CNA(i,j) = N_matrix(i,j)*(one - exp(-Z(i,j)))*(F(i,j)/Z(i,j));
    }
  } 
//compute beginning of year pop num@len (NL), catch@len (CL), and survey num@len (NLs);
  NL = pla*N_matrix;        
  // for CL should use pla for mid year;
  CL = pla*CNA;
  // for Nls should use pla for time of year, sf;
  NLs = pla*Ns;
  matrix<Type> log_NLs = log(NLs.array());

//compute biomass, catch biomass and ssb @ len, and total biomass and ssb;
  matrix<Type> B_matrix = weight.array()*NL.array();         
  matrix<Type> CB_matrix = weight.array()*CL.array();  
  matrix<Type> SSB_matrix = mat.array()*B_matrix.array(); 
  biomass = B_matrix.colwise().sum();
  ssb = SSB_matrix.colwise().sum();   
  vector<Type> catch_biomass = CB_matrix.colwise().sum();  
  vector<Type> harvest_rate = catch_biomass/biomass;
  
  // ssb relative to average, just FYI;
  vector<Type> rssb = Y*ssb/sum(ssb);
  vector<Type> q = exp(logq);

  //calculate model predicted survey index and residuals
  for(i = 0;i < n;++i){ 
    E_index(i) = 0.0;
    for(j = FRV_ilen1(i);j <=FRV_ilen2(i);++j){E_index(i) += q(j)*NLs(j,FRV_iyear(i));}
   }
   Elog_index = log(E_index);
   resid_index = log_index - Elog_index;  
   std_resid_index = resid_index/std_index;

  //calculate model predicted catch and residuals  
  for(i = 0;i < nc;++i){ 
    E_catch(i) = 0.0;
    for(j = CL_ilen1(i);j <=CL_ilen2(i);++j){E_catch(i) += CL(j,CL_iyear(i));}
   }
   Elog_catch = log(E_catch);
   resid_catch = log_catch - Elog_catch;
   std_resid_catch = resid_catch/std_catch; 

  //negative log-likelihoods
  //index
 	nll -= dnorm(resid_index,zero,std_index,true).sum();
  
  //catch
 	nll -= dnorm(resid_catch,zero,std_catch,true).sum(); 

  //logq 
  for(i = 3;i < 45;++i){nll -= dnorm(logq(i),logq(i-1),std_logq,true);} 

  //recruitment deviations 
  nll += SCALE(AR1(phi_logR),std_log_R)(log_Rec_dev);
    
  //Log_F_dev nll
  
  array<Type> log_F_dev1(A-3,Y); // A-3 because F_1 and F_2 = 0 , and F_A = F_A-1;
  for(i = 0;i < A-3;++i){ 
    for(j = 0;j < Y;++j){log_F_dev1(i,j) = log_F_dev(i,j);}
  }   
  
  nll += SCALE(SEPARABLE(AR1(phi_F_year),AR1(phi_F_age)),std_log_F)(log_F_dev1);

  //process error negative loglikelihood term
  //  for(int j = 0;j < Y-1;++j){
  //    for(int i = 0;i < A-1;++i){  
  //      nll -= dnorm(pe(i,j),zero,std_pe,true);       
  //    }
  //  }   	
  
  REPORT(pla);
  REPORT(mean_len); 
  REPORT(NL); 
  REPORT(NLs);     
  REPORT(CNA);    
  REPORT(CL); 
  REPORT(N_matrix); 
  REPORT(Ns);  
  REPORT(B_matrix);
  REPORT(SSB_matrix);  
  REPORT(CB_matrix);
  REPORT(Z);          
  REPORT(F);
  REPORT(biomass); 
  REPORT(catch_biomass); 
  REPORT(harvest_rate);
  REPORT(ssb);    
  REPORT(rssb);       
  REPORT(Elog_index);  
  REPORT(E_index);   
  REPORT(resid_index);  
  REPORT(std_resid_index); 
  REPORT(E_catch);  
  REPORT(Elog_catch); 
  REPORT(resid_catch);  
  REPORT(std_resid_catch); 
  REPORT(phi_logR);
   
  REPORT(log_meanR);  
  REPORT(std_log_R);
  REPORT(std_index);  
  REPORT(std_log_F); 
  REPORT(std_pe);
  REPORT(phi_logR);  
  REPORT(phi_F_age);  
  REPORT(phi_F_year);  
  REPORT(Linf);
  REPORT(vbk);  
  REPORT(po);   
  REPORT(len_o);
  REPORT(cv_len); 
  REPORT(log_Rec_dev); 
  REPORT(log_Rec);
  REPORT(pe);
  REPORT(log_F_main);
  REPORT(log_F_dev); 
  REPORT(logq);      
  REPORT(std_len);  
  REPORT(mean_len); 
  
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
