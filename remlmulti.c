/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


double eigen_invert(double *mat, int length, double *mat2, double *mat3, int flag)
{
int i, j, lwork, info;
double alpha, beta, wkopt, *work, det=0;

lwork=-1;
dsyev_("V", "U", &length, mat, &length, mat2, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &length, mat, &length, mat2, work, &lwork, &info);
if(info!=0){printf("Decomp error 2\n\n");exit(1);}
free(work);

if(mat3!=NULL)
{
for(j=0;j<length;j++)
{
if(mat2[j]>0)
{
det+=log(mat2[j]);
for(i=0;i<length;i++){mat3[i+j*length]=mat[i+j*length]*pow(mat2[j],-.5);}
}
else
{
det+=log(0.000001);
if(flag==1){printf("Decomp error 3\n");}
for(i=0;i<length;i++){mat3[i+j*length]=0;}
}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &length, &alpha, mat3, &length, mat3, &length, &beta, mat, &length);
}

return(det);
}	//end of eigen_invert

///////////////////////////

int multi_reml(double *Y, double *Z, double *X, int ns, int num_covars, int num_regs, int g, int *Xstarts, int *Xends, int *Xrec, int *Xrev, double *Xsums, int num_kins, double **mkins, int shortcut, double *U, double *E, double *reml, char *outfile, char **ids1, char **ids2, int rlength, double *rcentres, double *rmults, float *rweights, char **rprednames, char *ral1, char *ral2, double *gmults, float *gweights, char **gprednames, float adjust, float adjust2, double likenull, double *Mnull)
{
//g=-1 is reml, g=0 is null, else g>0 for genes
int i, i2, j, k, k2, r, count;
int nfree, one=1, nlost, token;
double alpha, beta, sum, value;

int gflag, wnum, *herfixes;
double gam, gam2a, gam2b, gam3, like, like2, like3, d2, dd2;
double *AI, *AI2, *AI3, *AIsave, *BI, *J, *G;
double *lambdas, *lamdiffs, *deltas, *hers, statlrt, pvalrt, statscore, pvascore;

double detV, *ZTVZ, *ZTVZ2, *ZTVZ3, detZTVZ, *PY, *PPY, **KPY, **PKPY, *traces;
double *ZTVY, *Yadj, *VYadj, *thetas, **mGtemp, **mG, **mG2, *effects, **effectsfull;
double *kintemp, *kintemp2, *UTG, *UTG2, *XTX, *XTX2, *XTX3, *XTG;
double *V, *V2, *V3, *VZ, *VZZTVZ, *P, *PX, *XTPY;
double *UTX, *DUTX, *XTVX, *XTVX2, *XTVX3, detXTVX, *DUTXXTVX;
double *UTY, *UTZ, *D, detD, detC, *F, *FUTZ, *FUTZZTVZ, *H, *HUTY, *HHUTY, *UTKPY, *HUTKPY;	//???
double *UTYadj, *FUTYadj, *XTVYadj;

FILE *output1, *output2, *output3, *output4, *output5;
char filename1[500], filename2[500], filename3[500], filename4[500];


//set gflag and save gene X
gflag=0;if(g>0){gflag=1;}
if(gflag==1)	//fill reml assuming gene does not contribute, check gene length and save predictors
{
reml[0]=0;reml[1]=0;reml[2]=likenull;reml[3]=0;reml[4]=1;reml[5]=0;reml[6]=1;
reml[7]=log(1.0/99999);reml[8]=0;reml[9]=-10;
if(Xstarts[num_regs]==Xends[num_regs]){return(0);}

token=Xends[num_regs]-Xstarts[num_regs];
G=malloc(sizeof(double)*ns*token);
for(i=0;i<ns*token;i++){G[i]=X[Xstarts[num_regs]*ns+i];}
}

nfree=ns-num_covars;
wnum=0;if(num_regs+gflag>0){wnum=Xends[num_regs+gflag-1];}

if(num_regs+gflag>0)	//sort things related to X
{
if(adjust>0&&adjust2>=adjust)
{printf("Error, the total weight for kinship 1 (%.2f) should be larger than the sum of total weights for each regions (%.2f)\n", adjust, adjust2);exit(1);}

if(shortcut==1)	//allocate and set things associated with X
{
UTX=malloc(sizeof(double)*ns*wnum);
DUTX=malloc(sizeof(double)*ns*wnum);
XTVX=malloc(sizeof(double)*wnum*wnum);
XTVX2=malloc(sizeof(double)*wnum);
XTVX3=malloc(sizeof(double)*wnum*wnum);
DUTXXTVX=malloc(sizeof(double)*ns*wnum);

if(num_kins==0)
{
for(i=0;i<ns*wnum;i++){UTX[i]=X[i];}
}
else
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &ns, &wnum, &ns, &alpha, U, &ns, X, &ns, &beta, UTX, &ns);
}
}
}	//end of sorting X

////////

ZTVZ=malloc(sizeof(double)*num_covars*num_covars);
ZTVZ2=malloc(sizeof(double)*num_covars);
ZTVZ3=malloc(sizeof(double)*num_covars*num_covars);
PY=malloc(sizeof(double)*ns);
PPY=malloc(sizeof(double)*ns);

KPY=malloc(sizeof(double*)*(num_kins+num_regs+gflag));
PKPY=malloc(sizeof(double*)*(num_kins+num_regs+gflag));
for(k=0;k<num_kins+num_regs+gflag;k++)
{KPY[k]=malloc(sizeof(double)*ns);
PKPY[k]=malloc(sizeof(double)*ns);}

traces=malloc(sizeof(double)*(num_kins+num_regs+gflag));
AI=malloc(sizeof(double)*(num_kins+num_regs+gflag)*(num_kins+num_regs+gflag));
AI2=malloc(sizeof(double)*(num_kins+num_regs+gflag));
AI3=malloc(sizeof(double)*(num_kins+num_regs+gflag)*(num_kins+num_regs+gflag));
AIsave=malloc(sizeof(double)*(num_kins+num_regs+gflag)*(num_kins+num_regs+gflag));
BI=malloc(sizeof(double)*(num_kins+num_regs+gflag));
J=malloc(sizeof(double)*(num_kins+num_regs+gflag)*(num_kins+num_regs+gflag));

ZTVY=malloc(sizeof(double)*num_covars);
Yadj=malloc(sizeof(double)*ns);
VYadj=malloc(sizeof(double)*ns);
thetas=malloc(sizeof(double)*num_covars);

if(num_kins+num_regs+gflag>0)
{
mG=malloc(sizeof(double*)*(num_kins+num_regs+gflag));
mGtemp=malloc(sizeof(double*)*(num_kins+num_regs+gflag));
mG2=malloc(sizeof(double*)*(num_kins+num_regs+gflag));
for(k=0;k<num_kins+num_regs+gflag;k++)
{
mG[k]=malloc(sizeof(double)*ns);
mGtemp[k]=malloc(sizeof(double)*ns);
mG2[k]=malloc(sizeof(double)*ns);	//might not use all of these
}
}

if(shortcut==0)	
{
V=malloc(sizeof(double)*ns*ns);
V2=malloc(sizeof(double)*ns);
V3=malloc(sizeof(double)*ns*ns);
VZ=malloc(sizeof(double)*num_covars*ns);
VZZTVZ=malloc(sizeof(double)*ns*num_covars);
P=malloc(sizeof(double)*ns*ns);

if(num_kins>0)
{kintemp=malloc(sizeof(double)*ns*ns);kintemp2=malloc(sizeof(double)*ns);}
}

if(shortcut==1)
{
UTY=malloc(sizeof(double)*ns);
UTZ=malloc(sizeof(double)*ns*num_covars);
if(num_kins==0)
{
for(i=0;i<ns;i++)
{
UTY[i]=Y[i];
for(j=0;j<num_covars;j++){UTZ[i+j*ns]=Z[i+j*ns];}
}
}
else
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &ns, &num_covars, &ns, &alpha, U, &ns, Z, &ns, &beta, UTZ, &ns);
dgemv_("T", &ns, &ns, &alpha, U, &ns, Y, &one, &beta, UTY, &one);
}

D=malloc(sizeof(double)*ns);
F=malloc(sizeof(double)*ns*ns);
FUTZ=malloc(sizeof(double)*ns*num_covars);
FUTZZTVZ=malloc(sizeof(double)*ns*ns);
H=malloc(sizeof(double)*ns*ns);
HUTY=malloc(sizeof(double)*ns);
HHUTY=malloc(sizeof(double)*ns);
UTKPY=malloc(sizeof(double)*ns);
HUTKPY=malloc(sizeof(double)*ns);

UTYadj=malloc(sizeof(double)*ns);
FUTYadj=malloc(sizeof(double)*ns);
}

///////////////////////////

//allocate her related variables
lambdas=malloc(sizeof(double)*(num_kins+num_regs+gflag));
lamdiffs=malloc(sizeof(double)*(num_kins+num_regs+gflag));
deltas=malloc(sizeof(double)*(num_kins+num_regs+gflag));
hers=malloc(sizeof(double)*(num_kins+num_regs+gflag));
herfixes=malloc(sizeof(int)*(num_kins+num_regs+gflag));

//set starting heritabilities - perhaps use priors ???
for(k=0;k<num_kins+num_regs;k++){hers[k]=.5/(num_kins+num_regs);}
if(gflag==1){hers[num_kins+num_regs]=.01;}
for(k=0;k<num_kins+num_regs+gflag;k++){herfixes[k]=0;}
nlost=0;

//turn these into lambdas then deltas - adjust>0 => num_kins=1
sum=0;for(k=0;k<num_kins+num_regs+gflag;k++){sum+=hers[k];}
for(k=0;k<num_kins+num_regs+gflag;k++)
{
lambdas[k]=hers[k]/(1-sum);
if(adjust==0){deltas[k]=lambdas[k];}
if(adjust>0&&k==0){deltas[0]=adjust/(adjust-adjust2)*lambdas[0];}
if(adjust>0&&k>0){deltas[k]=lambdas[k]-Xsums[k-1]/(adjust-adjust2)*lambdas[0];}
}

///////////////////////////

//now solve
#include "remlsolve.c"

///////////////////////////

if(num_kins+num_regs+gflag>0)	//for sds get jacobian - dlambdas/dhers - row k is dls/dhk
{
token=num_kins+num_regs+gflag;
sum=0;for(k=0;k<token;k++){sum+=hers[k];}
for(k=0;k<num_kins+num_regs+gflag;k++)
{
for(k2=0;k2<token;k2++){J[k+k2*token]=hers[k2]*pow(1-sum,-2);}
J[k+k*token]+=pow(1-sum,-1);
}

//want AI = J * AIsave * t(J), with zeros added, then invert
alpha=1.0;beta=0.0;
dgemm_("N", "T", &token, &token, &token, &alpha, AIsave, &token, J, &token, &beta, AI3, &token);
dgemm_("N", "N", &token, &token, &token, &alpha, J, &token, AI3, &token, &beta, AI, &token);
for(k=0;k<token;k++)
{
if(herfixes[k]>=2)
{
for(k2=0;k2<token;k2++){AI[k+k2*token]=0;AI[k2+k*token]=0;}
AI[k+k*token]=1;
}
}
(void)eigen_invert(AI, num_kins+num_regs+gflag, AI2, AI3, 0);
for(k=0;k<token;k++)
{
if(herfixes[k]>=2){AI[k+k*(num_kins+num_regs+gflag)]=0;}
}
}

/////////

//get fixed effects, print, then breeding values and region effects

//start with ZTinvVY
alpha=1.0;beta=0.0;
if(shortcut==0)
{dgemv_("T", &ns, &num_covars, &alpha, VZ, &ns, Y, &one, &beta, ZTVY, &one);}
else
{dgemv_("T", &ns, &num_covars, &alpha, FUTZ, &ns, UTY, &one, &beta, ZTVY, &one);}

//fixed effects are (ZTinvVZ)^-1 t(invVZ) Y with variance matrix (ZTinvVZ)^-1
alpha=1.0;beta=0.0;
dgemv_("T", &num_covars, &num_covars, &alpha, ZTVZ, &num_covars, ZTVY, &one, &beta, thetas, &one);

if(g==-1)	//print results
{
sprintf(filename2,"%s.reml", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fprintf(output2, "Kins:%d Regions:%d\n", num_kins, num_regs);
if(num_kins+num_regs>0){fprintf(output2, "Blupfile: %s.indi.blp\n", outfile);}
else{fprintf(output2, "Blupfile: none\n");}
if(num_regs>0){fprintf(output2, "Regfile: %s.reg.blup\n", outfile);}
else{fprintf(output2, "Regfile: none\n");}
//??? did have a line for whether completed properly or not
fprintf(output2, "Component Heritability SD\n");
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f\n", k+1, hers[k], pow(AI[k+k*(num_kins+num_regs+gflag)],.5));}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f\n", r+1, hers[num_kins+r], pow(AI[num_kins+r+(num_kins+r)*(num_kins+num_regs+gflag)],.5));}
fprintf(output2, "Intercept %.4f %.4f\n", thetas[0], pow(ZTVZ[0],.5));
for(j=1;j<num_covars;j++){fprintf(output2, "Fixed%d %.6f %.6f\n",j, thetas[j], pow(ZTVZ[j+j*num_covars],.5));} 
fclose(output2);
}

//breeding values are kin_k lambda_k invV (Y-fixed) = kin_k lambda_k UFUT (Y-fixed)

for(i=0;i<ns;i++){Yadj[i]=Y[i];}
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_covars, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

alpha=1.0;beta=0.0;
if(shortcut==0)
{
dgemv_("N", &ns, &ns, &alpha, V, &ns, Yadj, &one, &beta, VYadj, &one);
}
else
{
if(num_kins==0)
{
dgemv_("N", &ns, &ns, &alpha, F, &ns, Yadj, &one, &beta, VYadj, &one);
}
else
{
dgemv_("T", &ns, &ns, &alpha, U, &ns, Yadj, &one, &beta, UTYadj, &one);
dgemv_("N", &ns, &ns, &alpha, F, &ns, UTYadj, &one, &beta, FUTYadj, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, FUTYadj, &one, &beta, VYadj, &one);
}
}

//get kin_k invV (Y-fixed) (save multiplying by lambda[k] until afterwards)
for(k=0;k<num_kins;k++)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, mkins[k], &ns, VYadj, &one, &beta, mGtemp[k], &one);
}
for(r=0;r<num_regs+gflag;r++)
{
token=Xends[r]-Xstarts[r];
XTVYadj=malloc(sizeof(double)*token);
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("T", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, VYadj, &one, &beta, XTVYadj, &one);
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, XTVYadj, &one, &beta, mGtemp[num_kins+r], &one);
free(XTVYadj);
}

if(adjust>0)	//mG[0] should become N/N*mG[0] - N1/N*mG[1] - ...
{
for(i=0;i<ns;i++)
{
mGtemp[0][i]=mGtemp[0][i]*adjust/(adjust-adjust2);
for(r=0;r<num_regs+gflag;r++){mGtemp[0][i]-=Xsums[r]/(adjust-adjust2)*mGtemp[1+r][i];}
}
}

//now multiply by lambdas
for(k=0;k<num_kins+num_regs+gflag;k++)
{
for(i=0;i<ns;i++)
{
if(hers[k]>0){mG[k][i]=lambdas[k]*mGtemp[k][i];}
else{mG[k][i]=0.0;}
}
}

/////////

//for kinships, also want to get inv(kin) mG = U E^-1 (UT mG)
for(k=0;k<num_kins;k++)
{
if(hers[k]==0)
{
for(i=0;i<ns;i++){mG2[k][i]=0;}
}
else
{
if(shortcut==0)		//then do them all manually
{
if(adjust==0)
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){kintemp[i+i2*ns]=mkins[k][i+i2*ns];}
}
}
else	//must scale and subtract regions (must have k=0)
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){kintemp[i+i2*ns]=adjust/(adjust-adjust2)*mkins[0][i+i2*ns];}
}
for(r=0;r<num_regs+gflag;r++)
{
token=Xends[r]-Xstarts[r];
alpha=-1.0/(adjust-adjust2);
beta=1.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, kintemp, &ns);
}
}

(void)eigen_invert(kintemp, ns, kintemp2, NULL, 0);
UTG=malloc(sizeof(double)*ns);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &ns, &alpha, kintemp, &ns, mG[k], &one, &beta, UTG, &one);
for(i=0;i<ns;i++)
{
if(kintemp2[i]>0){UTG[i]=UTG[i]/kintemp2[i];}
else{UTG[i]=0;}
}
dgemv_("N", &ns, &ns, &alpha, kintemp, &ns, UTG, &one, &beta, mG2[k], &one);
free(UTG);
}

if(shortcut==1)	//then must have k=0
{
UTG=malloc(sizeof(double)*ns);
UTG2=malloc(sizeof(double)*ns);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &ns, &alpha, U, &ns, mG[k], &one, &beta, UTG, &one);

if(adjust==0)	//UTG2 = invE UTG
{
for(i=0;i<ns;i++)
{
if(E[i]>0){UTG2[i]=UTG[i]/E[i];}
else{UTG2[i]=0;}
}
}
else	//UTG2 = D' - D'UTX inv(-W1I+XTUD'UTX) XTUD' where D'=W1/WinvE
{
for(i=0;i<ns;i++)
{
if(E[i]>0){D[i]=(adjust-adjust2)/adjust/E[i];}
else{D[i]=0;}
}

//get XTVX = XTU invD' UTX - W1
for(i=0;i<ns;i++)
{
for(j=0;j<wnum;j++){DUTX[i+j*ns]=UTX[i+j*ns]*D[i];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &wnum, &wnum, &ns, &alpha, UTX, &ns, DUTX, &ns, &beta, XTVX, &wnum);
for(j=0;j<wnum;j++){XTVX[j+j*wnum]-=(adjust-adjust2);}

(void)eigen_invert(XTVX, wnum, XTVX2, XTVX3, 0);

//then pre multiply by invD'UTX and post by t(invD'XTU), subtracting from D'
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &wnum, &wnum, &alpha, DUTX, &ns, XTVX, &wnum, &beta, DUTXXTVX, &ns);
alpha=-1.0;beta=0.0;
dgemm_("N", "T", &ns, &ns, &wnum, &alpha, DUTXXTVX, &ns, DUTX, &ns, &beta, F, &ns);
for(i=0;i<ns;i++){F[i+i*ns]+=D[i];}

alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, F, &ns, UTG, &one, &beta, UTG2, &one);
}	//end of adjust>0

alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, U, &ns, UTG2, &one, &beta, mG2[k], &one);
free(UTG);free(UTG2);
}	//end of shortcut=1
}}	//end of her[k]>0 and k loop

/////////

if(g==-1)	//save
{
sprintf(filename3,"%s.indi.blp", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error opening %s\n\n",filename3);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output3, "%s\t%s\t", ids1[i], ids2[i]);
for(k=0;k<num_kins;k++){fprintf(output3, "%f\t%f\t", mG2[k][i], mG[k][i]);}
for(r=0;r<num_regs;r++){fprintf(output3, "0\t%f\t", mG[num_kins+r][i]);}
fprintf(output3, "\n");
}
fclose(output3);
}

///////////////////////////

if(num_regs+gflag>0)	//get region (and gene) effect sizes and save
{
effects=malloc(sizeof(double)*wnum);
for(j=0;j<wnum;j++){effects[j]=0;}
if(num_regs>0)
{
effectsfull=malloc(sizeof(double*)*(num_regs+1));
for(r=0;r<num_regs+1;r++)
{
effectsfull[r]=malloc(sizeof(double)*rlength);
for(j=0;j<rlength;j++){effectsfull[r][j]=0;}
}
}

for(r=0;r<num_regs+gflag;r++)
{
if(hers[num_kins+r]>0)	//when her=0 leave at zero
{
token=Xends[r]-Xstarts[r];
if(token<=ns)	//then effects = inv(XTX) XTg
{
XTX=malloc(sizeof(double)*token*token);
XTX2=malloc(sizeof(double)*token);
XTX3=malloc(sizeof(double)*token*token);
XTG=malloc(sizeof(double)*token);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &token, &token, &ns, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, XTX, &token);
(void)eigen_invert(XTX, token, XTX2, XTX3, 0);
dgemv_("T", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, mG[num_kins+r], &one, &beta, XTG, &one);
dgemv_("N", &token, &token, &alpha, XTX, &token, XTG, &one, &beta, effects+Xstarts[r], &one);
free(XTX);free(XTX2);free(XTX3);free(XTG);
}	//end of token<=ns

if(token>ns)	//then effects = XT inv(XXT) g
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, kintemp, &ns);
(void)eigen_invert(kintemp, ns, kintemp2, NULL, 0);
UTG=malloc(sizeof(double)*ns);
dgemv_("T", &ns, &ns, &alpha, kintemp, &ns, mG[num_kins+r], &one, &beta, UTG, &one);
for(i=0;i<ns;i++)
{
if(kintemp2[i]>0){UTG[i]=UTG[i]/kintemp2[i];}
else{UTG[i]=0;}
}
dgemv_("N", &ns, &ns, &alpha, kintemp, &ns, UTG, &one, &beta, mG2[num_kins+r], &one);
dgemv_("T", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, mG2[num_kins+r], &one, &beta, effects+Xstarts[r], &one);
free(UTG);
}	//end of token>ns

if(r<num_regs)
{
for(j=Xstarts[r];j<Xends[r];j++)
{effectsfull[r][Xrev[j]]+=effects[j]*rmults[Xrev[j]]*pow(rweights[Xrev[j]],.5);
effectsfull[num_regs][Xrev[j]]++;}
}
}}	//end of her>0 and r loop

if(g==-1)	//print these
{
sprintf(filename4,"%s.reg.blup", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error opening %s\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tA1\tA2\tCentre\t");
for(r=0;r<num_regs;r++){fprintf(output4, "Region%d\t", r+1);}
fprintf(output4, "\n");
for(j=0;j<rlength;j++)
{
if(effectsfull[num_regs][j]>0)
{
fprintf(output4,"%s\t%c\t%c\t%f\t", rprednames[j], ral1[j], ral2[j], rcentres[j]);
for(r=0;r<num_regs;r++){fprintf(output4, "%f\t", effectsfull[r][j]);}
fprintf(output4,"\n");
}
}
fclose(output4);
}

if(g>0)	//print gene effects
{
if((output5=fopen(outfile,"a"))==NULL)
{printf("Error opening %s\n\n",outfile);exit(1);}
for(j=Xstarts[num_regs];j<Xends[num_regs];j++)
{
fprintf(output5,"%d %s %f\n", g, gprednames[Xrev[j]], effects[j]*gmults[Xrev[j]]*pow(gweights[Xrev[j]],.5));
}
fclose(output5);
}

free(effects);
if(num_regs>0){for(r=0;r<num_regs+1;r++){free(effectsfull[r]);}free(effectsfull);}
}	//end of num_regs>0 & g>0

/////////

if(g==0)	//need to save P or H and like
{
if(shortcut==0)
{
for(i=0;i<ns*ns;i++){Mnull[i]=P[i];}
}
if(shortcut==1)
{
for(i=0;i<ns*ns;i++){Mnull[i]=H[i];}
}
reml[2]=like;
}

if(g>0)	//do lrt and score tests and load up reml
{
statlrt=2*(like-likenull);
pvalrt=cdfN(-pow(statlrt,.5));
if(statlrt<0){pvalrt=1;}

//get bf ???

token=Xends[num_regs]-Xstarts[num_regs];
if(adjust==0)	//use nullP or nullH and get PY
{
if(shortcut==0)
{
for(i=0;i<ns*ns;i++){P[i]=Mnull[i];}
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, Y, &one, &beta, PY, &one);
}
if(shortcut==1)
{
for(i=0;i<ns*ns;i++){H[i]=Mnull[i];}
if(num_kins==0)	//then U is the identity
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, Y, &one, &beta, PY, &one);
}
else
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, UTY, &one, &beta, HUTY, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, HUTY, &one, &beta, PY, &one);
}
}
}	//end of adjust=0
else	//need to reiterate with her set to zero
{
if(herfixes[num_kins+num_regs]<2)	//re-iterate with gene her zero
{
lambdas[num_kins+num_regs]=0;hers[num_kins+num_regs]=0;nlost++;
for(j=Xstarts[num_regs];j<Xends[num_regs];j++)
{
for(i=0;i<ns;i++)
{
X[i+j*ns]=0;
if(shortcut==1){UTX[i+j*ns]=0;}
}}
#include "remlsolve.c"
}
}

if(herfixes[num_kins+num_regs]>=2)	//restore X and perhaps UTX
{
for(i=0;i<ns*token;i++){X[Xstarts[num_regs]*ns+i]=G[i];}
if(shortcut==1)
{
if(num_kins==0)
{
for(i=0;i<ns*token;i++){UTX[Xstarts[num_regs]*ns+i]=X[Xstarts[num_regs]*ns+i];}
}
else
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &ns, &token, &ns, &alpha, U, &ns, X+Xstarts[num_regs]*ns, &ns, &beta, UTX+Xstarts[num_regs]*ns, &ns);
}
}
}

//get KPY
XTPY=malloc(sizeof(double)*token);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &token, &alpha, X+Xstarts[num_regs]*ns, &ns, PY, &one, &beta, XTPY, &one);
alpha=1.0/Xsums[num_regs];beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+Xstarts[num_regs]*ns, &ns, XTPY, &one, &beta, KPY[num_kins+num_regs], &one);
free(XTPY);

//and PKPY
if(shortcut==0)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, KPY[num_kins+num_regs], &one, &beta, PKPY[num_kins+num_regs], &one);
}
if(shortcut==1)	//U H UT KPY
{
alpha=1.0;beta=0.0;
if(num_kins==0)	//then U is diagonal
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, KPY[num_kins+num_regs], &one, &beta, PKPY[num_kins+num_regs], &one);
}
else
{
dgemv_("T", &ns, &ns, &alpha, U, &ns, KPY[num_kins+num_regs], &one, &beta, UTKPY, &one);
dgemv_("N", &ns, &ns, &alpha, H, &ns, UTKPY, &one, &beta, HUTKPY, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, HUTKPY, &one, &beta, PKPY[num_kins+num_regs], &one);
}
}

//get tr(PK)
traces[num_kins+num_regs]=0;
PX=malloc(sizeof(double)*ns*token);
if(shortcut==0)	//tr(PK) = tr(PX XT/N)
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &token, &ns, &alpha, P, &ns, X+Xstarts[num_regs]*ns, &ns, &beta, PX, &ns);
for(i=0;i<ns;i++)
{
for(j=0;j<token;j++)
{traces[num_kins+num_regs]+=PX[i+j*ns]*X[i+(Xstarts[num_regs]+j)*ns]/Xsums[num_regs];}
}
}
if(shortcut==1)	//tr(PK) = tr(HUTX XTU) (can store HUTX in PX)
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &token, &ns, &alpha, H, &ns, UTX+Xstarts[num_regs]*ns, &ns, &beta, PX, &ns);
for(i=0;i<ns;i++)
{
for(j=0;j<token;j++)
{traces[num_kins+num_regs]+=PX[i+j*ns]*UTX[i+(Xstarts[num_regs]+j)*ns]/Xsums[num_regs];}
}
}
free(PX);

gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
gam2a=0;for(i=0;i<ns;i++){gam2a+=PKPY[num_kins+num_regs][i]*Y[i];}
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[num_kins+num_regs][i]*KPY[num_kins+num_regs][i];}

d2=.5*nfree/gam*gam2a-.5*traces[num_kins+num_regs];
dd2=-.5*nfree/gam*(gam3-pow(gam2a,2)/gam);
statscore=d2/pow(-dd2,.5);
pvascore=cdfN(-statscore);
if(dd2>0){pvascore=1.0;}

reml[0]=hers[num_kins+num_regs];
reml[1]=pow(AI[num_kins+num_regs+(num_kins+num_regs)*(num_kins+num_regs+gflag)],.5);
reml[2]=like;reml[3]=statlrt;reml[4]=pvalrt/2;reml[5]=statscore;reml[6]=pvascore;
//??? do something for log and bf
}	//end of g>0

/////////////////////////// 

//free variables

free(lambdas);free(lamdiffs);free(deltas);free(hers);free(herfixes);

if(num_regs+gflag>0)
{
free(G);
if(shortcut==1){free(UTX);free(DUTX);free(XTVX);free(XTVX2);free(XTVX3);free(DUTXXTVX);}
}

free(ZTVZ);free(ZTVZ2);free(ZTVZ3);free(PY);free(PPY);free(traces);
for(k=0;k<num_kins+num_regs+gflag;k++){free(KPY[k]);free(PKPY[k]);}
if(num_kins+num_regs+gflag>0){free(KPY);free(PKPY);}
free(AI);free(AI2);free(AI3);free(AIsave);free(BI);free(J);
free(ZTVY);free(Yadj);free(VYadj);free(thetas);

for(k=0;k<num_kins+num_regs+gflag;k++){free(mG[k]);free(mGtemp[k]);free(mG2[k]);}
if(num_kins+num_regs+gflag>0){free(mG);free(mGtemp);free(mG2);}

if(shortcut==0)
{
free(V);free(V2);free(V3);free(VZ);free(VZZTVZ);free(P);
if(num_kins>0){free(kintemp);free(kintemp2);}
}

if(shortcut==1)
{
free(UTY);free(UTZ);free(D);free(F);free(FUTZ);free(FUTZZTVZ);
free(H);free(HUTY);free(HHUTY);free(UTKPY);free(HUTKPY);
free(UTYadj);free(FUTYadj);
}

return(0);
}	//end of multi_reml





/*
X=malloc(sizeof(double)*ns*(maxlength+glength));
Xsums=malloc(sizeof(double)*(num_regs+gflag));
Xstarts=malloc(sizeof(int)*(num_regs+gflag));
Xends=malloc(sizeof(int)*(num_regs+gflag));
Xrec=malloc(sizeof(int)*(maxlength+glength));
Xrev=malloc(sizeof(int)*(maxlength+glength));

wnum=0;adjust2=0;
for(r=0;r<num_regs;r++)
{
Xstarts[r]=wnum;Xsums[r]=0;
for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
if(rmults[j2]!=-1&&rweights[j2]>0)
{
for(i=0;i<ns;i++){X[i+wnum*ns]=rdata[kindex[1+i]+j2*ns];}
Xrec[wnum]=r;
Xrev[wnum]=j2;
Xsums[r]+=rweights[j2];
wnum++;
}}
Xends[r]=wnum;
adjust2+=Xsums[r];
}
if(gflag==1)
{
Xstarts[num_regs]=wnum;Xsums[num_regs]=0;
for(j=0;j<glength;j++)
{
if(gmults[j]!=-1&&gweights[j]>0)
{
for(i=0;i<ns;i++){X[i+wnum*ns]=gdata[kindex[1+i]+j*ns];}
Xrec[wnum]=num_regs;
Xrev[wnum]=j;
Xsums[num_regs]+=gweights[j];
wnum++;
}}
Xends[num_regs]=wnum;
adjust2+=Xsums[num_regs];
}
*/
