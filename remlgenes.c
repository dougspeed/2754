/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

int prepare_gene_reml(double *Y, double *Z, double *YTY, double *ZTY, double *ZTZ, double *detZTZ, double *YTCY, double *resp, int *kindex, int num_samples_use, double *covar, int num_covars)
{
int i, j, k;
int ns, one=1, info;
double alpha, beta;

ns=kindex[0];

//fill Y
for(i=0;i<ns;i++){Y[i]=resp[kindex[1+i]];}

//calc YTY
*YTY=0;for(i=0;i<ns;i++){*YTY+=pow(Y[i],2);}

//fill Z
for(i=0;i<ns;i++)
{
for(j=0;j<num_covars;j++){Z[i+j*ns]=covar[kindex[1+i]+j*num_samples_use];}
}

//calc ZTY, ZTZ, detZTZ, invZTZ (stored in ZTZ)
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_covars, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, covar, &ns, covar, &ns, &beta, ZTZ, &num_covars);

dpotrf_("U", &num_covars, ZTZ, &num_covars, &info);
if(info!=0)
{printf("Error with covariates 1\n\n");exit(1);}
*detZTZ=0;for(j=0;j<num_covars;j++){*detZTZ+=2*log(ZTZ[j+j*num_covars]);}
dpotri_("U", &num_covars, ZTZ, &num_covars, &info);
if(info!=0)
{printf("Error with covariates 2\n\n");exit(1);}

for(j=0;j<num_covars;j++)
{
for(k=0;k<j;k++){ZTZ[j+k*num_covars]=ZTZ[k+j*num_covars];}
}

//get YTCY = YTY - YTZ (inv)ZTZ ZTY and null likelihood
*YTCY=*YTY;
for(j=0;j<num_covars;j++)
{
for(k=0;k<num_covars;k++){*YTCY-=ZTY[j]*ZTZ[j+k*num_covars]*ZTY[k];}
}

return(0);
}	//end of prepare_gene_reml


///////////////////////////

int gene_reml(double *reml, double *Y, double *Z, double YTY, double *ZTY, double *ZTZ, double detZTZ, double YTCY, int num_covars, double *data, int *kindex, int num_samples_use, int length, double *mults, float *weights, float priora, float priorb)
{
int i, j, k, count;
int ns, nfree, one=1, info, lwork, best, stop;
double alpha, beta, wkopt, *work;

int wnum;
double wsum, lbetaab, S1, S2, S3, T1, T2, T3;
double gam, like, like2, like3, maxlike, likenull, deriv, dderiv, d2, dd2;
double lambda, lamdiff, her, hernew, sd, sdlog;
double statlrt, pvalrt, statscore, pvascore, remlbf;
double *X, *XTY, *XTZ, *XTCX, *XTCXtemp, *E, *U, *D, *Dtemp1, *Dtemp2;


ns=kindex[0];
nfree=ns-num_covars;
likenull=-.5*nfree*(1+log(2*M_PI*YTCY/nfree))-.5*detZTZ;

//allocations
X=malloc(sizeof(double)*ns*length);

//fill X
wnum=0;wsum=0;
for(j=0;j<length;j++)
{
if(mults[j]!=-1&&weights[j]>0)
{
for(i=0;i<ns;i++){X[i+wnum*ns]=data[kindex[1+i]+j*num_samples_use];}
wnum++;wsum+=weights[j];
}	//end of use this predictor
}

//fill reml assuming gene not tested
reml[0]=(double)wnum;reml[1]=wsum;reml[2]=0;reml[3]=0;
reml[4]=0;reml[5]=1;reml[6]=0;reml[7]=1;
reml[8]=log(1.0/99999);reml[9]=0;reml[10]=-10;reml[11]=YTCY/nfree;

if(wsum>0)
{
XTY=malloc(sizeof(double)*wnum);
XTZ=malloc(sizeof(double)*wnum*num_covars);
XTCX=malloc(sizeof(double)*wnum*wnum);
XTCXtemp=malloc(sizeof(double)*wnum*num_covars);
E=malloc(sizeof(double)*wnum);
U=malloc(sizeof(double)*wnum*wnum);
D=malloc(sizeof(double)*wnum);

//set lbetaab (use 1/wsum if not using priors)
lbetaab=-log(wsum);
if(priora!=-1){lbetaab=lgamma(priora)+lgamma(priorb)-lgamma(priora+priorb);}

//calc XTY, XTZ, XTCX = XTX - XTZ (inv)ZTZ * ZTX
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &wnum, &alpha, X, &ns, Y, &one, &beta, XTY, &one);
dgemm_("T", "N", &wnum, &num_covars, &ns, &alpha, X, &ns, Z, &ns, &beta, XTZ, &wnum);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &wnum, &wnum, &ns, &alpha, X, &ns, X, &ns, &beta, XTCX, &wnum);
dgemm_("N", "N", &wnum, &num_covars, &num_covars, &alpha, XTZ, &wnum, ZTZ, &num_covars, &beta, XTCXtemp, &wnum);
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &wnum, &wnum, &num_covars, &alpha, XTCXtemp, &wnum, XTZ, &wnum, &beta, XTCX, &wnum);

//decomp XTCX
for(j=0;j<wnum;j++)
{
for(k=0;k<wnum;k++){U[j+k*wnum]=XTCX[j+k*wnum];}
}

lwork=-1;
dsyev_("V", "U", &wnum, U, &wnum, E, &wkopt, &lwork, &info );
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &wnum, U, &wnum, E, work, &lwork, &info);
free(work);
if(info!=0){printf("error\n");}

//get D=UTXTCY = UT XTY - UT XTZ (inv)ZTZ ZTY
Dtemp1=malloc(sizeof(double)*num_covars);
Dtemp2=malloc(sizeof(double)*wnum);

alpha=1.0;beta=0.0;
dgemv_("T", &wnum, &wnum, &alpha, U, &wnum, XTY, &one, &beta, D, &one);
dgemv_("N", &num_covars, &num_covars, &alpha, ZTZ, &num_covars, ZTY, &one, &beta, Dtemp1, &one);
dgemv_("N", &wnum, &num_covars, &alpha, XTZ, &wnum, Dtemp1, &one, &beta, Dtemp2, &one);

alpha=-1.0;beta=1.0;
dgemv_("T", &wnum, &wnum, &alpha, U, &wnum, Dtemp2, &one, &beta, D, &one);

free(Dtemp1);
free(Dtemp2);


//ready for REML - adding on prior alog(w) +(b-1)log(lam) -(a+b)log(w+lam) -log(beta) 
//first test lambdas .1W, 1W, 10W, 100W, 1000W
best=-9;
for(k=-2;k<4;k++)
{
lambda=wsum*pow(10,k);
S1=0;for(j=0;j<wnum;j++){if(E[j]+lambda>0){S1+=log(E[j]+lambda);}}
T1=0;for(j=0;j<wnum;j++){T1+=pow(D[j],2)*pow(E[j]+lambda,-1);}
gam=YTCY-T1;
like=-.5*nfree*(1+log(2*M_PI*gam/nfree))-.5*S1+.5*wnum*log(lambda)-.5*detZTZ;
like+=priora*log(wsum)+(priorb-1)*log(lambda)-(priora+priorb)*log(wsum+lambda)-lbetaab;

if(best==-9){best=k;maxlike=like;}
if(like>maxlike){best=k;maxlike=like;}
}

lambda=wsum*pow(10,best);
her=wsum/(wsum+lambda);

like2=0;stop=0;count=0;
while(1)
{
S1=0;for(j=0;j<wnum;j++){if(E[j]+lambda>0){S1+=log(E[j]+lambda);}}
S2=0;for(j=0;j<wnum;j++){S2+=pow(E[j]+lambda,-1);}
S3=0;for(j=0;j<wnum;j++){S3+=pow(E[j]+lambda,-2);}
T1=0;for(j=0;j<wnum;j++){T1+=pow(D[j],2)*pow(E[j]+lambda,-1);}
T2=0;for(j=0;j<wnum;j++){T2+=pow(D[j],2)*pow(E[j]+lambda,-2);}
T3=0;for(j=0;j<wnum;j++){T3+=pow(D[j],2)*pow(E[j]+lambda,-3);}

//get gamma then derivs and like
gam=YTCY-T1;
deriv=-.5*nfree/gam*T2-.5*S2+.5*wnum/lambda +(priorb-1)/lambda-(priora+priorb)/(wsum+lambda);
dderiv=.5*nfree/gam*(pow(T2,2)/gam+2*T3)+.5*S3-.5*wnum*pow(lambda,-2) -(priorb-1)*pow(lambda,-2)+(priora+priorb)*pow(wsum+lambda,-2);
like3=like2;like2=like;
like=-.5*nfree*(1+log(2*M_PI*gam/nfree))-.5*S1+.5*wnum*log(lambda)-.5*detZTZ;
like+=priora*log(wsum)+(priorb-1)*log(lambda)-(priora+priorb)*log(wsum+lambda)-lbetaab;

//always want to break before updating
if(stop==1){break;}
if(count==1000){printf("Gene did not finish after %d REML iterations, I hope it got close\n", count);break;}

lamdiff=deriv/dderiv;
if(lamdiff>lambda-0.00001*wsum){lamdiff=lambda-0.00001*wsum;}	//this implies h near 1
if(lamdiff<lambda-99999*wsum){lamdiff=lambda-99999*wsum;}	//implies h near 0

lambda=lambda-lamdiff;
hernew=wsum/(wsum+lambda);
if(hernew-her>0.1){hernew=her+.1;lambda=wsum*(1-hernew)/hernew;}
if(her-hernew>0.1){hernew=her-.1;lambda=wsum*(1-hernew)/hernew;}
her=hernew;

if(abs(like-like2)<0.0001&&abs(like2-like3)<0.0001){stop=1;}
count++;
}

if(lambda==0.00001*wsum){her=1.0;lambda=0;}
if(lambda==99999*wsum){her=0;like=likenull;}

//get sd and mltest
dd2=pow(wsum+lambda,4)*pow(wsum,-2)*dderiv;
sd=pow(-dd2,-.5);
if(dd2>0){sd=-1;}
if(her==0){sd=0;}

statlrt=2*(like-likenull);
pvalrt=cdfN(-pow(statlrt,.5));
if(statlrt<0){pvalrt=1;}

//get bf from integral
remlbf=(like+.5*log(-2*M_PI/dd2)-likenull)/log(10);
if(dd2>0){remlbf=-10;}

//for score test get deriv for 1/lambda at zero
S2=0;for(j=0;j<wnum;j++){S2+=E[j];}
S3=0;for(j=0;j<wnum;j++){S3+=pow(E[j],2);}
T2=0;for(j=0;j<wnum;j++){T2+=pow(D[j],2);}
T3=0;for(j=0;j<wnum;j++){T3+=pow(D[j],2)*E[j];}

d2=.5*nfree*T2/YTCY-.5*S2;
dd2=.5*nfree/YTCY*(pow(T2,2)/YTCY-2*T3)+.5*S3;
statscore=d2/pow(-dd2,.5);
if(statscore>0){statscore=-statscore;}
pvascore=2*cdfN(statscore);
if(dd2>0){pvascore=1.0;}
	
//want also sd for log (1/lambda) = log (h/(1-h)/w)
dd2=pow(lambda,2)*dderiv;
sdlog=pow(-dd2,-.5);

reml[0]=(double)wnum;reml[1]=wsum;reml[2]=her;reml[3]=sd;
reml[4]=statlrt;reml[5]=pvalrt/2;reml[6]=statscore;reml[7]=pvascore;
if(her>0){reml[8]=log(wsum)-log(lambda);reml[9]=sdlog;reml[10]=remlbf;reml[11]=gam/nfree;}

free(XTY);free(XTZ);free(XTCX);free(XTCXtemp);free(E);free(U);free(D);
}	//end of if wnum>0

free(X);

return(0);
}	//end of gene_reml




