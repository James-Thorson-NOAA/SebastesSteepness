DATA_SECTION
  init_int m;  //# Number of stocks

  init_int knots;                 //# Number of evaluations per stock
  init_vector h_profX(1,knots);      //# Location of evaluations per stock

  //#THIS IS NEW CODE
  init_vector StockWeight(1,m);
  //# END NEW CODE

  init_matrix h_profY(1,m,1,knots);   //# NLL for each stock and evaluation

  //#THIS IS NEW CODE
  init_int TestVal;
  !! if (TestVal != 123456) { cout << "Test Number is not 123456" << endl; exit(1); }
  !! if(TestVal == 123456) { cout << "Data loaded" << endl; }
  //# END NEW CODE
  
  int mccounter;
  !! mccounter = 1;

PARAMETER_SECTION
  init_bounded_number mu_h(0.20001,0.99999);          //# Mean
  init_bounded_number tau(0.001,10);       //# LogSD

  vector h(1,m);
  vector h_prof(1,m);                 //# 
  vector beta(1,knots);               //# Logit-transformed steepness
  vector Jacob(1,knots);              //# Jacobian
  number mu; 

  matrix LikeMat(1,m,1,knots);      //# Record likelihood of each h value 

  number nu
  number s2

  //number v;
  //number s;
  //number delta;

  objective_function_value obj_fun;
  sdreport_number log_tau;                  //# SD

PRELIMINARY_CALCS_SECTION
  int Istk,Ih;
  double maxLik;
  cout << setprecision(4);
  //v = 10.0;
  //s = sqrt(.5);

  //# informatative prior on square(tau)
  nu = 10.;
  s2 = 0.5;
  //# Vague prior
  //  nu = 2.;
  //  s2 = 2.;
  //# Vaguer prior
  //  nu = 1.;
  //  s2 = .5;

  //# Calculate Logit-transformed X
  for (Ih=1;Ih<=knots;Ih++){
    if (h_profX(Ih) < 1 & h_profX(Ih) > 0.2){
      beta(Ih) = log((h_profX(Ih)-0.2)/(1-h_profX(Ih)));  //# Logit-transformed X
    }
    if (h_profX(Ih)==1){  
      beta(Ih) = log((0.999-0.2)/(1-0.999));
    }
    if (h_profX(Ih)==0.2){ 
      beta(Ih) = log((0.201-0.2)/(1-0.201));
    }
  }
  //# Delta is used to ensure that the integral works correctly
   //delta = (h_profX(2)-h_profX(1))/2.506628;
  //# Jacobian (which results in an equivalent expression to standard logit-normal)
  for (Ih=1;Ih<=knots;Ih++)
   Jacob(Ih) = 0.8*exp(beta(Ih))/square(1.0+exp(beta(Ih)));   //# reciprocal of derivative of logit-function = derivative of logistic function 
   
  //# Convert to rescaled Likelihood
  for (Istk=1;Istk<=m;Istk++)
   {
    maxLik = h_profY(Istk,1);
    for (Ih=2;Ih<=knots;Ih++)
     if (h_profY(Istk,Ih) < maxLik) maxLik = h_profY(Istk,Ih);
    for (Ih=1;Ih<=knots;Ih++)
     h_profY(Istk,Ih) = mfexp(maxLik-h_profY(Istk,Ih));
   }

PROCEDURE_SECTION
  int Istk, Ih;
  dvariable LikeStock,likeV;
  dvariable pi;

  //#THIS IS NEW CODE
  mu = log((mu_h-0.2)/0.8) - log(1-(mu_h-0.2)/0.8);
  //# END NEW CODE
  pi = 3.141592;
  obj_fun = 0;
  log_tau = log(tau);
  
  //# Summation across profile likelihoods for MLE (i.e. requires Jacobian) 
  for (Istk=1;Istk<=m;Istk++)
   {
    LikeStock = 0;
    for (Ih=1;Ih<=knots;Ih++){              
     LikeMat(Istk,Ih) = h_profY(Istk,Ih) * (1/sqrt(2*pi*square(tau)))*mfexp(-1*(square(beta(Ih)-mu)/(2*square(tau)))) / Jacob(Ih);
     LikeStock = LikeStock + LikeMat(Istk,Ih);
    }
    obj_fun -= StockWeight(Istk) * log(LikeStock); 
   }
 
  //#THIS IS NEW CODE
  //# Inverse chi-square prior on square(tau) 
  obj_fun += (((nu/2)+1)*log(square(tau)))+((nu*s2)/(2.0*square(tau)));    //# negative log-likelihood of 2nd definition of inverse-chi-squared from wikipedia with constant removed
  //# END NEW CODE

  //# Track 'Empirical Bayes' estimates of h and report parameters 
  if (mceval_phase()) 
   {
    for (Istk=1;Istk<=m;Istk++)
     {
      LikeStock = 0;
      for (Ih=1;Ih<=knots;Ih++)
       {
        likeV = h_profY(Istk,Ih)*exp(-0.5*square(beta(Ih)-mu)/square(tau))/tau;
        if (likeV > LikeStock) { h(Istk) = h_profX(Ih); LikeStock = likeV; }
       }
     }  
   
    ofstream hout("h.out", ios::app);
    hout << h << endl;
    hout.close();
    ofstream tauout("tau.out", ios::app);
    tauout << tau << endl;
    tauout.close();
    ofstream muout("mu.out", ios::app);
    muout << mu << endl;
    muout.close();
    ofstream h_profout("h_prof.out", ios::app);
    h_profout << h_prof << endl;
    h_profout.close();
   }

REPORT_SECTION
  report << "# mu\n" << mu << endl;
  report << "# tau\n" << tau << endl;
  report << "# m\n" << m << endl;
  report << "# knots\n" << knots << endl;
  report << "# beta\n" << beta << endl;
  report << "# Jacob\n" << Jacob << endl;
  report << "# h_profX\n" << h_profX << endl;
  report << "# h_profY\n" << h_profY << endl;
  report << "# LikeMat\n" << LikeMat << endl;
