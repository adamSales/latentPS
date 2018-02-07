data{
//Sample sizes
 int<lower=1> ncov;
 int<lower=1> nstudTO;
 int<lower=1> nstudTM;
 int<lower=2> nstudC;
 int<lower=1> nteacher;
 int<lower=1> nschool;
 int<lower=1> npair;

// indices
 int<lower=1,upper=nteacher> teacherTO[nstudTO];
 int<lower=1,upper=npair> pairTO[nstudTO];
 int<lower=1,upper=nschool> schoolTO[nstudTO];
 int<lower=1,upper=nteacher> teacherTM[nstudTM];
 int<lower=1,upper=npair> pairTM[nstudTM];
 int<lower=1,upper=nschool> schoolTM[nstudTM];
 int<lower=1,upper=nteacher> teacherC[nstudC];
 int<lower=1,upper=npair> pairC[nstudC];
 int<lower=1,upper=nschool> schoolC[nstudC];

// data data
 matrix[nstudTO,ncov] XtO;
 matrix[nstudTM,ncov] XtM;
 matrix[nstudC,ncov] Xc;
 real YtO[nstudTO];
 real YtM[nstudTM];
 real Yc[nstudC];
 vector<lower=0,upper=1>[nstudTO] MbarTO;

}
parameters{

 real alphaU;
 vector[ncov] betaU;
 vector[ncov] betaY;

 real a0;
 real a1;
 real b0;
 real b1;

 real teacherEffY[nteacher];
 real teacherEffU[nteacher];
 real pairEffect[npair];
 real schoolEffY[nschool];

 real<lower=0> sigTchY;
 real<lower=0> sigSclY;
 real<lower=0> sigY[2];
 real<lower=0> sigTchU;
 real<lower=0> sigU;

 vector<lower=0,upper=1>[nstudC] MbarC;
 vector<lower=0,upper=1>[nstudTM] MbarTM;
}

model{
 vector[nstudTO] muYtO;
 vector[nstudTM] muYtM;
 vector[nstudC] muYc;
 vector[nstudTO] muUtO;
 vector[nstudTM] muUtM;
 vector[nstudC] muUc;

 for(i in 1:nstudTO){
  muYtO[i]=teacherEffY[teacherTO[i]]+schoolEffY[schoolTO[i]]+pairEffect[pairTO[i]]+a0+b0+(a1+b1)*MbarTO[i];
  muUtO[i]=teacherEffU[teacherTO[i]];
 }

 for(i in 1:nstudTM){
  muYtM[i]=teacherEffY[teacherTM[i]]+schoolEffY[schoolTM[i]]+pairEffect[pairTM[i]]+a0+b0+(a1+b1)*MbarTM[i];
  muUtM[i]=teacherEffU[teacherTM[i]];
 }


 for(i in 1:nstudC){
  muYc[i]=teacherEffY[teacherC[i]]+schoolEffY[schoolC[i]]+pairEffect[pairC[i]]+a0+a1*MbarC[i];
  muUc[i]=teacherEffU[teacherC[i]];
 }

 //priors
 betaY~normal(0,2);
 betaU~normal(0,2);
 pairEffect~normal(0,2);

 a0~normal(0,1);
 a1~normal(0,1);
 b0~normal(0,1);
 b1~normal(0,1);

 // actual model:
 schoolEffY~normal(0,sigSclY);
 teacherEffU~normal(0,sigTchU);
 teacherEffY~normal(0,sigTchY);

 MbarTO~normal(alphaU+muUtO+XtO*betaU,sigU);
 MbarTM~normal(alphaU+muUtM+XtM*betaU,sigU);
 MbarC~normal(alphaU+muUc+Xc*betaU,sigU);
 YtO~normal(muYtO+XtO*betaY,sigY[2]);
 YtM~normal(muYtM+XtM*betaY,sigY[2]);
 Yc~normal(muYc+Xc*betaY,sigY[1]);
}