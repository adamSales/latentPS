data{
//Sample sizes
 int<lower=1> nsecWorked;
 int<lower=1> ncov;
 int<lower=1> nstud;
 int<lower=1> nteacher;
 int<lower=1> nsec;
 int<lower=1> nschool;
 int<lower=1> npair;

// indices
 int<lower=1,upper=nteacher> teacher[nstud];
 int<lower=1,upper=npair> pair[nstud];
 int<lower=1,upper=nschool> school[nstud];
 int<lower=1,upper=nstud> studentM[nsecWorked];
 int<lower=1,upper=nsec> section[nsecWorked];

// data data
 int<lower=0,upper=1> grad[nsecWorked];
 matrix[nstud,ncov] X;
 int<lower=0,upper=1> Z[nstud];
 real Y[nstud];

}
parameters{

 vector[nstud] studEff;

 vector[ncov] betaU;
 vector[ncov] betaY;

// real a0;
 real a1;
 real b0;
 real b1;

 real teacherEffY[nteacher];
// real teacherEffU[nteacher];
 real pairEffect[npair];
// real schoolEffU[nschool];
 real schoolEffY[nschool];
 real secEff[nsec];

 real<lower=0> sigTchY;
 real<lower=0> sigSclY;
 real<lower=0> sigY[2];
 //real<lower=0> sigTchU;
 //real<lower=0> sigSclU;
 real<lower=0> sigU;
}

model{
 real linPred[nsecWorked];
 vector[nstud] muY;
 //vector[nstud] muU;
 real useEff[nstud];
 real trtEff[nstud];
 real sigYI[nstud];


// grad model
 for(i in 1:nsecWorked)
  linPred[i]= secEff[section[i]]+studEff[studentM[i]];

 for(i in 1:nstud){
  useEff[i]=a1*studEff[i];
  trtEff[i]=b0+b1*studEff[i];
  //muU[i]=teacherEffU[teacher[i]]+schoolEffU[school[i]];
  muY[i]=teacherEffY[teacher[i]]+schoolEffY[school[i]]+pairEffect[pair[i]]+useEff[i]+Z[i]*trtEff[i];
  sigYI[i]=sigY[Z[i]+1];
 }

 //priors
 betaY~normal(0,2);
 betaU~normal(0,2);
 pairEffect~normal(0,2);

// a0~normal(0,1);
 a1~normal(0,1);
 b0~normal(0,1);
 b1~normal(0,1);


 schoolEffY~normal(0,sigSclY);
 //schoolEffU~normal(0,sigSclU);
 //teacherEffU~normal(0,sigTchU);
 teacherEffY~normal(0,sigTchY);

 grad~bernoulli_logit(linPred);

 studEff~normal(X*betaU,sigU);
 Y~normal(muY+X*betaY,sigYI);
}