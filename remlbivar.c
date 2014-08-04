/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


void stack(double *mat, double a, double b, double c, int length)
{
int j, j2;
int token=2*length;

for(j=0;j<token;j++)
{
for(j2=0;j2<token;j2++){mat[j+j2*token]=0;}
}

for(j=0;j<length;j++)
{
mat[j+j*token]=a;
mat[length+j+j*token]=b;
mat[j+(length+j)*token]=b;
mat[length+j+(length+j)*token]=c;
}
}

///////

double inner(double *A, double *B, int ll)
{
int j;
double ret=0;
for(j=0;j<ll;j++){ret+=A[j]*B[j];}
return(ret);
}

double trace1(double *A, int ll)
{
int j;
double ret=0;
for(j=0;j<ll;j++){ret+=A[j+j*ll];}
return(ret);
}

double trace2(double *A, double *B, int ll)
{
int j, j2;
double ret=0;
for(j=0;j<ll;j++)
{
for(j2=0;j2<ll;j2++){ret+=A[j+j2*ll]*B[j2+j*ll];}
}
return(ret);
}

///////

int gene_bivar(double *bivar, double *Ya, double *Yb, double *Za, double *Zb, double *ZTYa, double *ZTYb, double *ZTZa, double *ZTZb, double detZTZa, double detZTZb, double YTCYa, double YTCYb, int num_covars, double *data, int **respindex, int num_samples_use, int length, double *mults, float *weights)
{
int i, j, j2, count;
int nsa, nsb, nfreea, nfreeb, token, one=1, five=5, info, lwork;
double alpha, beta, wkopt, *work;

int wnum;
double wsum, ea, eb, la, lb, ro, ronew;
double adds[5], diff, like, likenull, statlrt, pvalrt, statscore, pvascore;
double hera, herb, heranew, herbnew, herremla, sdremla=-1, herremlb, sdremlb=-1, dd2;

double *Xa, *Xb, *XTYa, *XTYb, *XTZa, *XTZb, *XTCXa, *XTCXb, *XTCXtemp;
double *Da, *Db, *Dtemp, *D, *T, detT, *T2, *T3, *TD;
double *Tea, *Teb, *Tla, *Tlb, *Tro, *TeaTD, *TebTD, *TlaTD, *TlbTD, *TroTD;
double *TTeaTD, *TTebTD, *TTlaTD, *TTlbTD, *TTroTD, *TTea, *TTeb, *TTla, *TTlb, *TTro;
double *Teaea, *Teaeb, *Teala, *Tealb, *Tearo, *TeaeaTD, *TeaebTD, *TealaTD, *TealbTD, *TearoTD;
double *Tebea, *Tebeb, *Tebla, *Teblb, *Tebro, *TebeaTD, *TebebTD, *TeblaTD, *TeblbTD, *TebroTD;
double *Tlaea, *Tlaeb, *Tlala, *Tlalb, *Tlaro, *TlaeaTD, *TlaebTD, *TlalaTD, *TlalbTD, *TlaroTD;
double *Tlbea, *Tlbeb, *Tlbla, *Tlblb, *Tlbro, *TlbeaTD, *TlbebTD, *TlblaTD, *TlblbTD, *TlbroTD;
double *Troea, *Troeb, *Trola, *Trolb, *Troro, *TroeaTD, *TroebTD, *TrolaTD, *TrolbTD, *TroroTD;

double dea, deb, dla, dlb, dro;
double deaea, deaeb, deala, dealb, dearo, debea, debeb, debla, deblb, debro;
double dlaea, dlaeb, dlala, dlalb, dlaro, dlbea, dlbeb, dlbla, dlblb, dlbro, droea, droeb, drola, drolb, droro;
double AI[25], AIsave[25], AI2[5], AI3[25], BI[5];


nsa=respindex[0][0];
nsb=respindex[1][0];
nfreea=nsa-num_covars;
nfreeb=nsb-num_covars;

//allocations
Xa=malloc(sizeof(double)*nsa*length);
Xb=malloc(sizeof(double)*nsb*length);
XTYa=malloc(sizeof(double)*length);
XTYb=malloc(sizeof(double)*length);
XTZa=malloc(sizeof(double)*length*num_covars);
XTZb=malloc(sizeof(double)*length*num_covars);
XTCXa=malloc(sizeof(double)*length*length);
XTCXb=malloc(sizeof(double)*length*length);
XTCXtemp=malloc(sizeof(double)*length*num_covars);

//fill Xs
wnum=0;wsum=0;
for(j=0;j<length;j++)
{
if(mults[j]!=-1&&weights[j]>0)
{
for(i=0;i<nsa;i++){Xa[i+wnum*nsa]=data[respindex[0][1+i]+j*num_samples_use];}
for(i=0;i<nsb;i++){Xb[i+wnum*nsb]=data[respindex[1][1+i]+j*num_samples_use];}
wnum++;wsum+=weights[j];
}	//end of use this predictor
}

//must have wsum>0 to be here

//calc XTYs, XTZs, XTCXs = XTXs - XTZs (inv)ZTZs * ZTXs
alpha=1.0;beta=0.0;
dgemv_("T", &nsa, &wnum, &alpha, Xa, &nsa, Ya, &one, &beta, XTYa, &one);
dgemm_("T", "N", &wnum, &num_covars, &nsa, &alpha, Xa, &nsa, Za, &nsa, &beta, XTZa, &wnum);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &wnum, &wnum, &nsa, &alpha, Xa, &nsa, Xa, &nsa, &beta, XTCXa, &wnum);
dgemm_("N", "N", &wnum, &num_covars, &num_covars, &alpha, XTZa, &wnum, ZTZa, &num_covars, &beta, XTCXtemp, &wnum);
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &wnum, &wnum, &num_covars, &alpha, XTCXtemp, &wnum, XTZa, &wnum, &beta, XTCXa, &wnum);

alpha=1.0;beta=0.0;
dgemv_("T", &nsb, &wnum, &alpha, Xb, &nsb, Yb, &one, &beta, XTYb, &one);
dgemm_("T", "N", &wnum, &num_covars, &nsb, &alpha, Xb, &nsb, Zb, &nsb, &beta, XTZb, &wnum);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &wnum, &wnum, &nsb, &alpha, Xb, &nsb, Xb, &nsb, &beta, XTCXb, &wnum);
dgemm_("N", "N", &wnum, &num_covars, &num_covars, &alpha, XTZb, &wnum, ZTZb, &num_covars, &beta, XTCXtemp, &wnum);
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &wnum, &wnum, &num_covars, &alpha, XTCXtemp, &wnum, XTZb, &wnum, &beta, XTCXb, &wnum);

free(XTCXtemp);

//get Ds=XTCYs, each of which = XTY - XTZ (inv)ZTZ ZTY, then load into D
Da=malloc(sizeof(double)*wnum);
Db=malloc(sizeof(double)*wnum);
Dtemp=malloc(sizeof(double)*num_covars);
D=malloc(sizeof(double)*2*wnum);

for(j=0;j<wnum;j++){Da[j]=XTYa[j];}
alpha=1.0;beta=0.0;
dgemv_("N", &num_covars, &num_covars, &alpha, ZTZa, &num_covars, ZTYa, &one, &beta, Dtemp, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &wnum, &num_covars, &alpha, XTZa, &wnum, Dtemp, &one, &beta, Da, &one);

for(j=0;j<wnum;j++){Db[j]=XTYb[j];}
alpha=1.0;beta=0.0;
dgemv_("N", &num_covars, &num_covars, &alpha, ZTZb, &num_covars, ZTYb, &one, &beta, Dtemp, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &wnum, &num_covars, &alpha, XTZb, &wnum, Dtemp, &one, &beta, Db, &one);

for(j=0;j<wnum;j++){D[j]=Da[j];D[wnum+j]=Db[j];}

free(Dtemp);

//set starting values
ea=bivar[12];la=wsum/bivar[13];eb=bivar[14];lb=wsum/bivar[15];ro=0.0;
hera=wsum/(wsum+ea*la);herb=wsum/(wsum+eb*lb);

//more allcoations
token=2*wnum;

T=malloc(sizeof(double)*token*token);
T2=malloc(sizeof(double)*token);
T3=malloc(sizeof(double)*token*token);
TD=malloc(sizeof(double)*token);

Tea=malloc(sizeof(double)*token*token);
Teb=malloc(sizeof(double)*token*token);
Tla=malloc(sizeof(double)*token*token);
Tlb=malloc(sizeof(double)*token*token);
Tro=malloc(sizeof(double)*token*token);

TeaTD=malloc(sizeof(double)*token);
TebTD=malloc(sizeof(double)*token);
TlaTD=malloc(sizeof(double)*token);
TlbTD=malloc(sizeof(double)*token);
TroTD=malloc(sizeof(double)*token);

TTeaTD=malloc(sizeof(double)*token);
TTebTD=malloc(sizeof(double)*token);
TTlaTD=malloc(sizeof(double)*token);
TTlbTD=malloc(sizeof(double)*token);
TTroTD=malloc(sizeof(double)*token);

TTea=malloc(sizeof(double)*token*token);
TTeb=malloc(sizeof(double)*token*token);
TTla=malloc(sizeof(double)*token*token);
TTlb=malloc(sizeof(double)*token*token);
TTro=malloc(sizeof(double)*token*token);

Teaea=malloc(sizeof(double)*token*token);
Teaeb=malloc(sizeof(double)*token*token);
Teala=malloc(sizeof(double)*token*token);
Tealb=malloc(sizeof(double)*token*token);
Tearo=malloc(sizeof(double)*token*token);

TeaeaTD=malloc(sizeof(double)*token*token);
TeaebTD=malloc(sizeof(double)*token*token);
TealaTD=malloc(sizeof(double)*token*token);
TealbTD=malloc(sizeof(double)*token*token);
TearoTD=malloc(sizeof(double)*token*token);

Tebea=malloc(sizeof(double)*token*token);
Tebeb=malloc(sizeof(double)*token*token);
Tebla=malloc(sizeof(double)*token*token);
Teblb=malloc(sizeof(double)*token*token);
Tebro=malloc(sizeof(double)*token*token);

TebeaTD=malloc(sizeof(double)*token*token);
TebebTD=malloc(sizeof(double)*token*token);
TeblaTD=malloc(sizeof(double)*token*token);
TeblbTD=malloc(sizeof(double)*token*token);
TebroTD=malloc(sizeof(double)*token*token);

Tlaea=malloc(sizeof(double)*token*token);
Tlaeb=malloc(sizeof(double)*token*token);
Tlala=malloc(sizeof(double)*token*token);
Tlalb=malloc(sizeof(double)*token*token);
Tlaro=malloc(sizeof(double)*token*token);

TlaeaTD=malloc(sizeof(double)*token*token);
TlaebTD=malloc(sizeof(double)*token*token);
TlalaTD=malloc(sizeof(double)*token*token);
TlalbTD=malloc(sizeof(double)*token*token);
TlaroTD=malloc(sizeof(double)*token*token);

Tlbea=malloc(sizeof(double)*token*token);
Tlbeb=malloc(sizeof(double)*token*token);
Tlbla=malloc(sizeof(double)*token*token);
Tlblb=malloc(sizeof(double)*token*token);
Tlbro=malloc(sizeof(double)*token*token);

TlbeaTD=malloc(sizeof(double)*token*token);
TlbebTD=malloc(sizeof(double)*token*token);
TlblaTD=malloc(sizeof(double)*token*token);
TlblbTD=malloc(sizeof(double)*token*token);
TlbroTD=malloc(sizeof(double)*token*token);

Troea=malloc(sizeof(double)*token*token);
Troeb=malloc(sizeof(double)*token*token);
Trola=malloc(sizeof(double)*token*token);
Trolb=malloc(sizeof(double)*token*token);
Troro=malloc(sizeof(double)*token*token);

TroeaTD=malloc(sizeof(double)*token*token);
TroebTD=malloc(sizeof(double)*token*token);
TrolaTD=malloc(sizeof(double)*token*token);
TrolbTD=malloc(sizeof(double)*token*token);
TroroTD=malloc(sizeof(double)*token*token);


diff=1;count=0;
while(1)
{
//now will want to store Tinv in T
for(j=0;j<wnum;j++)
{
for(j2=0;j2<wnum;j2++)
{
T[j+j2*token]=XTCXa[j+j2*wnum]*ea;
T[wnum+j+(wnum+j2)*token]=XTCXb[j+j2*wnum]*eb;
T[wnum+j+j2*token]=0;
T[j+(wnum+j2)*token]=0;
}
T[j+j*token]+=ea*ea*la;
T[wnum+j+(wnum+j)*token]+=eb*eb*lb;
T[wnum+j+j*token]=-ro*pow(la*lb,.5)*ea*eb;
T[j+(wnum+j)*token]=-ro*pow(la*lb,.5)*ea*eb;
}

//get Tinv and detT then Tinv D and like
lwork=-1;
dsyev_("V", "U", &token, T, &token, T2, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &token, T, &token, T2, work, &lwork, &info);
if(info!=0){printf("Decomp error 2\n\n");exit(1);}
free(work);

//set T3 to U E^-.5 (U currently stored in T, E in T2) 
detT=0;
for(j2=0;j2<token;j2++)
{
if(T2[j2]>0)
{
detT+=log(T2[j2]);
for(j=0;j<token;j++){T3[j+j2*token]=T[j+j2*token]*pow(T2[j2],-.5);}
}
else
{
detT+=log(.000001);printf("Error doug\n");
for(j=0;j<token;j++){T3[j+j2*token]=0;}
}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &token, &token, &token, &alpha, T3, &token, T3, &token, &beta, T, &token);
dgemv_("N", &token, &token, &alpha, T, &token, D, &one, &beta, TD, &one);

like=-.5*YTCYa/ea-.5*YTCYb/eb+(wnum-.5*nfreea)*log(ea)+(wnum-.5*nfreeb)*log(eb)-.5*detT+.5*wnum*log(la*lb*(1-ro*ro))-.5*(nfreea+nfreeb)*log(2*M_PI)-.5*detZTZa-.5*detZTZb;
for(j=0;j<token;j++){like+=.5*D[j]*TD[j];}

if(count==0){likenull=like;}

//get first derivs in full matrix form

stack(Tea, 2*ea*la, -ro*pow(la*lb,.5)*eb, 0, wnum);
stack(Teb, 0, -ro*pow(la*lb,.5)*ea, 2*eb*lb, wnum);
stack(Tla, ea*ea, -.5*ro*pow(lb/la,.5)*ea*eb, 0, wnum);
stack(Tlb, 0, -.5*ro*pow(la/lb,.5)*ea*eb, eb*eb, wnum);
stack(Tro, 0, -pow(la*lb,.5)*ea*eb, 0, wnum);

for(j=0;j<wnum;j++)
{
for(j2=0;j2<wnum;j2++)
{
Tea[j+j2*token]+=XTCXa[j+j2*wnum];
Teb[wnum+j+(wnum+j2)*token]+=XTCXb[j+j2*wnum];
}
}

//get first derivs * TinvD

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Tea, &token, TD, &one, &beta, TeaTD, &one);
dgemv_("N", &token, &token, &alpha, Teb, &token, TD, &one, &beta, TebTD, &one);
dgemv_("N", &token, &token, &alpha, Tla, &token, TD, &one, &beta, TlaTD, &one);
dgemv_("N", &token, &token, &alpha, Tlb, &token, TD, &one, &beta, TlbTD, &one);
dgemv_("N", &token, &token, &alpha, Tro, &token, TD, &one, &beta, TroTD, &one);

//and premult by Tinv

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, T, &token, TeaTD, &one, &beta, TTeaTD, &one);
dgemv_("N", &token, &token, &alpha, T, &token, TebTD, &one, &beta, TTebTD, &one);
dgemv_("N", &token, &token, &alpha, T, &token, TlaTD, &one, &beta, TTlaTD, &one);
dgemv_("N", &token, &token, &alpha, T, &token, TlbTD, &one, &beta, TTlbTD, &one);
dgemv_("N", &token, &token, &alpha, T, &token, TroTD, &one, &beta, TTroTD, &one);

//get Tinv * first derivs

alpha=1.0;beta=0.0;
dgemm_("N", "N", &token, &token, &token, &alpha, T, &token, Tea, &token, &beta, TTea, &token);
dgemm_("N", "N", &token, &token, &token, &alpha, T, &token, Teb, &token, &beta, TTeb, &token);
dgemm_("N", "N", &token, &token, &token, &alpha, T, &token, Tla, &token, &beta, TTla, &token);
dgemm_("N", "N", &token, &token, &token, &alpha, T, &token, Tlb, &token, &beta, TTlb, &token);
dgemm_("N", "N", &token, &token, &token, &alpha, T, &token, Tro, &token, &beta, TTro, &token);


//get second derivs in full matrix form and postmult by TinvD

stack(Teaea, 2*la, 0, 0, wnum);
stack(Teaeb, 0, -ro*pow(la*lb,.5), 0, wnum);
stack(Teala, 2*ea, -.5*ro*pow(lb/la,.5)*eb, 0, wnum);
stack(Tealb, 0, -.5*ro*pow(la/lb,.5)*eb, 0, wnum);
stack(Tearo, 0, -pow(la*lb,.5)*eb, 0, wnum);

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Teaea, &token, TD, &one, &beta, TeaeaTD, &one);
dgemv_("N", &token, &token, &alpha, Teaeb, &token, TD, &one, &beta, TeaebTD, &one);
dgemv_("N", &token, &token, &alpha, Teala, &token, TD, &one, &beta, TealaTD, &one);
dgemv_("N", &token, &token, &alpha, Tealb, &token, TD, &one, &beta, TealbTD, &one);
dgemv_("N", &token, &token, &alpha, Tearo, &token, TD, &one, &beta, TearoTD, &one);

stack(Tebea, 0, -ro*pow(la*lb,.5), 0, wnum);
stack(Tebeb, 0, 0, 2*lb, wnum);
stack(Tebla, 0, -.5*ro*pow(lb/la,.5)*ea, 0, wnum);
stack(Teblb, 0, -.5*ro*pow(la/lb,.5)*ea, 2*eb, wnum);
stack(Tebro, 0, -pow(la*lb,.5)*ea, 0, wnum);

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Tebea, &token, TD, &one, &beta, TebeaTD, &one);
dgemv_("N", &token, &token, &alpha, Tebeb, &token, TD, &one, &beta, TebebTD, &one);
dgemv_("N", &token, &token, &alpha, Tebla, &token, TD, &one, &beta, TeblaTD, &one);
dgemv_("N", &token, &token, &alpha, Teblb, &token, TD, &one, &beta, TeblbTD, &one);
dgemv_("N", &token, &token, &alpha, Tebro, &token, TD, &one, &beta, TebroTD, &one);

stack(Tlaea, 2*ea, -.5*ro*pow(lb/la,.5)*eb, 0, wnum);
stack(Tlaeb, 0, -.5*ro*pow(lb/la,.5)*ea, 0, wnum);
stack(Tlala, 0, .25*ro*pow(lb/la/la/la,.5)*ea*eb, 0, wnum);
stack(Tlalb, 0, -.25*ro*pow(1/lb/la,.5)*ea*eb, 0, wnum);
stack(Tlaro, 0, -.5*pow(lb/la,.5)*ea*eb, 0, wnum);

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Tlaea, &token, TD, &one, &beta, TlaeaTD, &one);
dgemv_("N", &token, &token, &alpha, Tlaeb, &token, TD, &one, &beta, TlaebTD, &one);
dgemv_("N", &token, &token, &alpha, Tlala, &token, TD, &one, &beta, TlalaTD, &one);
dgemv_("N", &token, &token, &alpha, Tlalb, &token, TD, &one, &beta, TlalbTD, &one);
dgemv_("N", &token, &token, &alpha, Tlaro, &token, TD, &one, &beta, TlaroTD, &one);

stack(Tlbea, 0, -.5*ro*pow(la/lb,.5)*eb, 0, wnum);
stack(Tlbeb, 0, -.5*ro*pow(la/lb,.5)*ea, 2*eb, wnum);
stack(Tlbla, 0, -.25*ro*pow(1/la/lb,.5)*ea*eb, 0, wnum);
stack(Tlblb, 0, .25*ro*pow(la/lb/lb/lb,.5)*ea*eb, 0, wnum);
stack(Tlbro, 0, -.5*pow(la/lb,.5)*ea*eb, 0, wnum);

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Tlbea, &token, TD, &one, &beta, TlbeaTD, &one);
dgemv_("N", &token, &token, &alpha, Tlbeb, &token, TD, &one, &beta, TlbebTD, &one);
dgemv_("N", &token, &token, &alpha, Tlbla, &token, TD, &one, &beta, TlblaTD, &one);
dgemv_("N", &token, &token, &alpha, Tlblb, &token, TD, &one, &beta, TlblbTD, &one);
dgemv_("N", &token, &token, &alpha, Tlbro, &token, TD, &one, &beta, TlbroTD, &one);

stack(Troea, 0, -pow(la*lb,.5)*eb, 0, wnum);
stack(Troeb, 0, -pow(la*lb,.5)*ea, 0, wnum);
stack(Trola, 0, -.5*pow(lb/la,.5)*ea*eb, 0, wnum);
stack(Trolb, 0, -.5*pow(la/lb,.5)*ea*eb, 0, wnum);
stack(Troro, 0, 0, 0, wnum);

alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, Troea, &token, TD, &one, &beta, TroeaTD, &one);
dgemv_("N", &token, &token, &alpha, Troeb, &token, TD, &one, &beta, TroebTD, &one);
dgemv_("N", &token, &token, &alpha, Trola, &token, TD, &one, &beta, TrolaTD, &one);
dgemv_("N", &token, &token, &alpha, Trolb, &token, TD, &one, &beta, TrolbTD, &one);
dgemv_("N", &token, &token, &alpha, Troro, &token, TD, &one, &beta, TroroTD, &one);

//put together twice first derivatives

dea=-inner(TeaTD,TD,token)-trace1(TTea,token)+YTCYa/ea/ea+(2.0*wnum-nfreea)/ea;
deb=-inner(TebTD,TD,token)-trace1(TTeb,token)+YTCYb/eb/eb+(2*wnum-nfreeb)/eb;
dla=-inner(TlaTD,TD,token)-trace1(TTla,token)+wnum/la;
dlb=-inner(TlbTD,TD,token)-trace1(TTlb,token)+wnum/lb;
dro=-inner(TroTD,TD,token)-trace1(TTro,token)-2*wnum*ro/(1-ro*ro);

BI[0]=dea;BI[1]=deb;BI[2]=dla;BI[3]=dlb;BI[4]=dro;

//now second derivatives

deaea=2*inner(TeaTD,TTeaTD,token)-inner(TeaeaTD,TD,token)+1*trace2(TTea,TTea,token)-trace2(Teaea,T,token) -2*YTCYa/ea/ea/ea-(2*wnum-nfreea)/ea/ea;
deaeb=2*inner(TeaTD,TTebTD,token)-inner(TeaebTD,TD,token)+1*trace2(TTea,TTeb,token)-trace2(Teaeb,T,token);
deala=2*inner(TeaTD,TTlaTD,token)-inner(TealaTD,TD,token)+1*trace2(TTea,TTla,token)-trace2(Teala,T,token);
dealb=2*inner(TeaTD,TTlbTD,token)-inner(TealbTD,TD,token)+1*trace2(TTea,TTlb,token)-trace2(Tealb,T,token);
dearo=2*inner(TeaTD,TTroTD,token)-inner(TearoTD,TD,token)+1*trace2(TTea,TTro,token)-trace2(Tearo,T,token);

debea=2*inner(TebTD,TTeaTD,token)-inner(TebeaTD,TD,token)+1*trace2(TTeb,TTea,token)-trace2(Tebea,T,token);
debeb=2*inner(TebTD,TTebTD,token)-inner(TebebTD,TD,token)+1*trace2(TTeb,TTeb,token)-trace2(Tebeb,T,token) -2*YTCYb/eb/eb/eb-(2*wnum-nfreeb)/eb/eb;
debla=2*inner(TebTD,TTlaTD,token)-inner(TeblaTD,TD,token)+1*trace2(TTeb,TTla,token)-trace2(Tebla,T,token);
deblb=2*inner(TebTD,TTlbTD,token)-inner(TeblbTD,TD,token)+1*trace2(TTeb,TTlb,token)-trace2(Teblb,T,token);
debro=2*inner(TebTD,TTroTD,token)-inner(TebroTD,TD,token)+1*trace2(TTeb,TTro,token)-trace2(Tebro,T,token);

dlaea=2*inner(TlaTD,TTeaTD,token)-inner(TlaeaTD,TD,token)+1*trace2(TTla,TTea,token)-trace2(Tlaea,T,token);
dlaeb=2*inner(TlaTD,TTebTD,token)-inner(TlaebTD,TD,token)+1*trace2(TTla,TTeb,token)-trace2(Tlaeb,T,token);
dlala=2*inner(TlaTD,TTlaTD,token)-inner(TlalaTD,TD,token)+1*trace2(TTla,TTla,token)-trace2(Tlala,T,token) -wnum/la/la;
dlalb=2*inner(TlaTD,TTlbTD,token)-inner(TlalbTD,TD,token)+1*trace2(TTla,TTlb,token)-trace2(Tlalb,T,token);
dlaro=2*inner(TlaTD,TTroTD,token)-inner(TlaroTD,TD,token)+1*trace2(TTla,TTro,token)-trace2(Tlaro,T,token);

dlbea=2*inner(TlbTD,TTeaTD,token)-inner(TlbeaTD,TD,token)+1*trace2(TTlb,TTea,token)-trace2(Tlbea,T,token);
dlbeb=2*inner(TlbTD,TTebTD,token)-inner(TlbebTD,TD,token)+1*trace2(TTlb,TTeb,token)-trace2(Tlbeb,T,token);
dlbla=2*inner(TlbTD,TTlaTD,token)-inner(TlblaTD,TD,token)+1*trace2(TTlb,TTla,token)-trace2(Tlbla,T,token);
dlblb=2*inner(TlbTD,TTlbTD,token)-inner(TlblbTD,TD,token)+1*trace2(TTlb,TTlb,token)-trace2(Tlblb,T,token) -wnum/lb/lb;
dlbro=2*inner(TlbTD,TTroTD,token)-inner(TlbroTD,TD,token)+1*trace2(TTlb,TTro,token)-trace2(Tlbro,T,token);

droea=2*inner(TroTD,TTeaTD,token)-inner(TroeaTD,TD,token)+1*trace2(TTro,TTea,token)-trace2(Troea,T,token);
droeb=2*inner(TroTD,TTebTD,token)-inner(TroebTD,TD,token)+1*trace2(TTro,TTeb,token)-trace2(Troeb,T,token);
drola=2*inner(TroTD,TTlaTD,token)-inner(TrolaTD,TD,token)+1*trace2(TTro,TTla,token)-trace2(Trola,T,token);
drolb=2*inner(TroTD,TTlbTD,token)-inner(TrolbTD,TD,token)+1*trace2(TTro,TTlb,token)-trace2(Trolb,T,token);
droro=2*inner(TroTD,TTroTD,token)-inner(TroroTD,TD,token)+1*trace2(TTro,TTro,token)-trace2(Troro,T,token) -2*wnum*(1+ro*ro)/pow(1-ro*ro,2);

AI[0]=-deaea;AI[1]=-deaeb;AI[2]=-deala;AI[3]=-dealb;AI[4]=-dearo;
AI[5]=-debea;AI[6]=-debeb;AI[7]=-debla;AI[8]=-deblb;AI[9]=-debro;
AI[10]=-dlaea;AI[11]=-dlaeb;AI[12]=-dlala;AI[13]=-dlalb;AI[14]=-dlaro;
AI[15]=-dlbea;AI[16]=-dlbeb;AI[17]=-dlbla;AI[18]=-dlblb;AI[19]=-dlbro;
AI[20]=-droea;AI[21]=-droeb;AI[22]=-drola;AI[23]=-drolb;AI[24]=-droro;

for(j=0;j<25;j++){AIsave[j]=AI[j];}

if(count==0)
{
statscore=BI[4]/pow(AI[24],.5);
if(statscore>0){statscore=-statscore;}
pvascore=2*cdfN(statscore);
if(AI[24]>0){statscore=0;}
}


//invert AI and multiply by BI
lwork=-1;
dsyev_("V", "U", &five, AI, &five, AI2, &wkopt, &lwork, &info );
if(info!=0){printf("Becomp error 1\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &five, AI, &five, AI2, work, &lwork, &info);
if(info!=0){printf("Becomp error 2\n\n");exit(1);}
free(work);

//set AI3 to U E^-.5
for(j2=0;j2<5;j2++)
{
for(j=0;j<5;j++)
{
if(AI2[j2]>0){AI3[j+j2*5]=AI[j+j2*5]*pow(AI2[j2],-.5);}
else{AI3[j+j2*5]=0;}
}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &five, &five, &five, &alpha, AI3, &five, AI3, &five, &beta, AI, &five);
dgemv_("N", &five, &five, &alpha, AI, &five, BI, &one, &beta, adds, &one);

//break before updating
if(diff<0.0001){break;}
if(count==500){printf("Bivariate did not finish after %d REML iterations, I hope it got close\n", count);break;}

//move parameters, limiting how much rho, ea, eb and hers can move

if(adds[0]>0.1*YTCYa/nfreea){adds[0]=0.1*YTCYa/nfreea;}
if(adds[0]<(-0.1)*YTCYa/nfreea){adds[0]=-0.1*YTCYa/nfreea;}
ea=ea+adds[0];if(ea<0.0001){ea=0.0001;}

if(adds[1]>0.1*YTCYb/nfreeb){adds[1]=0.1*YTCYb/nfreeb;}
if(adds[1]<-0.1*YTCYb/nfreeb){adds[1]=-0.1*YTCYb/nfreeb;}
eb=eb+adds[1];if(eb<0.0001){eb=0.0001;}

if(adds[4]>0.1){adds[4]=0.1;}if(adds[4]<-0.1){adds[4]=-0.1;}
ronew=ro+adds[4];if(ronew<-0.9999){ronew=-0.9999;}if(ronew>0.9999){ronew=0.9999;}

if(adds[2]>0.1*pow(wsum+la*(1-ro*ro)*ea,2)/(1-ro*ro)/ea/wsum)
{adds[2]=0.1*pow(wsum+la*(1-ro*ro)*ea,2)/(1-ro*ro)/ea/wsum;}
if(adds[2]<-0.1*pow(wsum+la*(1-ro*ro)*ea,2)/(1-ro*ro)/ea/wsum)
{adds[2]=-0.1*pow(wsum+la*(1-ro*ro)*ea,2)/(1-ro*ro)/ea/wsum;}
la=la+adds[2];if(la<0.0001){la=0.0001;}if(la>99999){la=99999;}
heranew=wsum/(wsum+ea*la*(1-ronew*ronew));

if(adds[3]>0.1*pow(wsum+lb*(1-ro*ro)*eb,2)/(1-ro*ro)/eb/wsum)
{adds[3]=0.1*pow(wsum+lb*(1-ro*ro)*eb,2)/(1-ro*ro)/eb/wsum;}
if(adds[3]<-0.1*pow(wsum+lb*(1-ro*ro)*eb,2)/(1-ro*ro)/eb/wsum)
{adds[3]=-0.1*pow(wsum+lb*(1-ro*ro)*eb,2)/(1-ro*ro)/eb/wsum;}
lb=lb+adds[3];if(lb<0.0001){lb=0.0001;}if(lb>99999){lb=99999;}
herbnew=wsum/(wsum+eb*lb*(1-ronew*ronew));

diff=fabs(heranew-hera)+fabs(herbnew-herb)+fabs(ronew-ro);
hera=heranew;herb=herbnew;ro=ronew;

count++;
}	//end of while loop

//have not worked out how to get sds

//get lrt
statlrt=2*(like-likenull);
pvalrt=cdfN(-pow(statlrt,.5));
if(statlrt<0){pvalrt=1;}

bivar[2]=ro;bivar[3]=pow(AI[24],.5);if(AI[24]<0){bivar[3]=-1;}
bivar[4]=herremla;bivar[5]=sdremla;bivar[6]=herremlb;bivar[7]=sdremlb;
bivar[8]=statlrt;bivar[9]=pvalrt;bivar[10]=statscore;bivar[11]=pvascore;

free(Xa);free(Xb);free(XTYa);free(XTYb);free(XTZa);free(XTZb);free(XTCXa);free(XTCXb);
free(Da);free(Db);free(D);free(T);free(T2);free(T3);free(TD);

free(Tea);free(Teb);free(Tla);free(Tlb);free(Tro);free(TeaTD);free(TebTD);free(TlaTD);free(TlbTD);free(TroTD);
free(TTeaTD);free(TTebTD);free(TTlaTD);free(TTlbTD);free(TTroTD);free(TTea);free(TTeb);free(TTla);free(TTlb);free(TTro);
free(Teaea);free(Teaeb);free(Teala);free(Tealb);free(Tearo);free(TeaeaTD);free(TeaebTD);free(TealaTD);free(TealbTD);free(TearoTD);
free(Tebea);free(Tebeb);free(Tebla);free(Teblb);free(Tebro);free(TebeaTD);free(TebebTD);free(TeblaTD);free(TeblbTD);free(TebroTD);
free(Tlaea);free(Tlaeb);free(Tlala);free(Tlalb);free(Tlaro);free(TlaeaTD);free(TlaebTD);free(TlalaTD);free(TlalbTD);free(TlaroTD);
free(Tlbea);free(Tlbeb);free(Tlbla);free(Tlblb);free(Tlbro);free(TlbeaTD);free(TlbebTD);free(TlblaTD);free(TlblbTD);free(TlbroTD);
free(Troea);free(Troeb);free(Trola);free(Trolb);free(Troro);free(TroeaTD);free(TroebTD);free(TrolaTD);free(TrolbTD);free(TroroTD);

return(0);
}	//end of gene_bivar





