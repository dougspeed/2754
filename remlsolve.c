/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


if(g==-1)	//will be screen and file printing progress
{
printf("Iter\t");
for(k=0;k<num_kins;k++){printf("Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){printf("Her_R%d\t", r+1);}
printf("Likelihood\tDifference\tTarget\tNo_Set_to_Zero\n");

sprintf(filename1,"%s.progress", outfile);
if((output1=fopen(filename1,"w"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
fprintf(output1, "Iter\t");
for(k=0;k<num_kins;k++){fprintf(output1, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output1, "Her_R%d\t", r+1);}
fprintf(output1, "Likelihood\tDifference\tTarget\tNo_Set_to_Zero\n");
fclose(output1);
}

like=1;like2=0;count=0;
while(1)
{
//first want invV = (V1 + V2)^-1; store in V or as U F UT

if(shortcut==0)	//get V=(I+delta_1K1 + delta_2 K2 + ...)^-1 directly
{
//start with diagonal matrix
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){V[i+i2*ns]=0.0;}
V[i+i*ns]=1.0;
}

//add on delta_k K for kinships
for(k=0;k<num_kins;k++)
{
if(deltas[k]!=0)
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){V[i+i2*ns]+=deltas[k]*mkins[k][i+i2*ns];}
}
}}

//and for regions and gene
for(r=0;r<num_regs+gflag;r++)
{
if(deltas[num_kins+r]!=0)
{
token=Xends[r]-Xstarts[r];
alpha=deltas[num_kins+r]/Xsums[r];
beta=1.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, V, &ns);
}
}

//invert V
detV=eigen_invert(V, ns, V2, V3, 1);
}	//end of shortcut=0

if(shortcut==1)	//get F=Dinv-Dinv UTX (invC+XTX)^-1 XTU invD)
{
//find D, where V1=UDUT - starts as I - entries should always be nonzero?
for(i=0;i<ns;i++){D[i]=1.0;}

if(num_kins==1)	//add on E delta_0
{
if(deltas[0]!=0)
{
for(i=0;i<ns;i++){D[i]+=deltas[0]*E[i];}
}
}
for(i=0;i<ns;i++){if(D[i]<=0){printf("warn %d %f\n", i+1, D[i]);}}

//for detV need detD, detC and detXTVX (last is actually det (invC + XTVX))
detD=0;for(i=0;i<ns;i++){detD+=log(D[i]);}
detC=0;detXTVX=0;

if(num_regs+gflag>0)	//now the regions - can replace zero deltas with anything nonzero
{
//get XTVX = XTU invDUTX + C^-1 where C = delta/W
for(i=0;i<ns;i++)
{
for(j=0;j<wnum;j++){DUTX[i+j*ns]=UTX[i+j*ns]/D[i];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &wnum, &wnum, &ns, &alpha, UTX, &ns, DUTX, &ns, &beta, XTVX, &wnum);

for(j=0;j<wnum;j++)
{
r=Xrec[j];
value=deltas[num_kins+r]/Xsums[r];
if(value==0){value=1.0;}
XTVX[j+j*wnum]+=1.0/value;
detC+=log(value);
}

//invert XTVX - something might go wrong when adjusting ???
detXTVX=eigen_invert(XTVX, wnum, XTVX2, XTVX3, 1);

//then pre multiply by invDUTX and post by t(invDXTU), subtracting from invD
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &wnum, &wnum, &alpha, DUTX, &ns, XTVX, &wnum, &beta, DUTXXTVX, &ns);
alpha=-1.0;beta=0.0;
dgemm_("N", "T", &ns, &ns, &wnum, &alpha, DUTXXTVX, &ns, DUTX, &ns, &beta, F, &ns);
for(i=0;i<ns;i++){F[i+i*ns]+=1.0/D[i];}
}	//end of have region+gene bit
else	//F is simply Dinv - add detV ???
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){F[i+i2*ns]=0;}
F[i+i*ns]=1.0/D[i];
}
}

detV=detD+detC+detXTVX;
}	//end of shortcut=1

///////////////////////////

//now want P = invV - invVZ invZTVZ ZTinvV or H = F - F UTZ invZTVZ ZTU F

if(shortcut==0)	//get invVZ then ZT invVZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_covars, &ns, &alpha, V, &ns, Z, &ns, &beta, VZ, &ns);
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, Z, &ns, VZ, &ns, &beta, ZTVZ, &num_covars);
}

if(shortcut==1)	//get FUTZ, then ZTU FUTZ (stored in ZTVZ) - could simplify if F diagonal
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_covars, &ns, &alpha, F, &ns, UTZ, &ns, &beta, FUTZ, &ns);
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, UTZ, &ns, FUTZ, &ns, &beta, ZTVZ, &num_covars);
}

//invert ZTinvVZ
detZTVZ=eigen_invert(ZTVZ, num_covars, ZTVZ2, ZTVZ3, 1);

if(shortcut==0)	//pre and post multiply by invVZ and ZTinvV, and subtract from invV
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){P[i+i2*ns]=V[i+i2*ns];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_covars, &num_covars, &alpha, VZ, &ns, ZTVZ, &num_covars, &beta, VZZTVZ, &ns);
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &ns, &ns, &num_covars, &alpha, VZZTVZ, &ns, VZ, &ns, &beta, P, &ns);
}

if(shortcut==1)	//pre and post multiply by FUTZ and ZTUF, and subtract from F
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){H[i+i2*ns]=F[i+i2*ns];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_covars, &num_covars, &alpha, FUTZ, &ns, ZTVZ, &num_covars, &beta, FUTZZTVZ, &ns);
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &ns, &ns, &num_covars, &alpha, FUTZZTVZ, &ns, FUTZ, &ns, &beta, H, &ns);
}

///////////////////////////

//get PY, PPY, gam and like

if(shortcut==0)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, Y, &one, &beta, PY, &one);
dgemv_("N", &ns, &ns, &alpha, P, &ns, PY, &one, &beta, PPY, &one);
}

if(shortcut==1)
{
alpha=1.0;beta=0.0;
if(num_kins==0)	//then U is identity
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, Y, &one, &beta, PY, &one);
dgemv_("N", &ns, &ns, &alpha, H, &ns, PY, &one, &beta, PPY, &one);
}
else
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, UTY, &one, &beta, HUTY, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, HUTY, &one, &beta, PY, &one);
dgemv_("N", &ns, &ns, &alpha, H, &ns, HUTY, &one, &beta, HHUTY, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, HHUTY, &one, &beta, PPY, &one);
}
}

//get gam=YTPY, from which can get likelihood
gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
like3=like2;like2=like;
like=-.5*nfree*(1+log(2*M_PI*gam/nfree))-.5*detV-.5*detZTVZ;

if(nlost==num_kins+num_regs+gflag&&g>=0){break;}
if(nlost==num_kins+num_regs+gflag){printf("All heritabiliities are zero\n");break;}
if(fabs(like-like2)<0.0001&&fabs(like2-like3)<0.0001){break;}
if(count==50){printf("Did not finish after %d REML iterations, but generally current values will be usable\n\n", count);break;}

if(g==-1)	//screen / file print progress
{
if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
for(k=0;k<num_kins;k++){printf("%.6f\t", hers[k]);}
for(r=0;r<num_regs;r++){printf("%.6f\t", hers[num_kins+r]);}
printf("%f\t", like);
if(count==0){printf("n/a\t\t%f\t0\n", 0.0001);}
else{printf("%f\t%f\t%d\n", like-like2, 0.0001, nlost);}

if((output1=fopen(filename1,"a"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
if(count==0){fprintf(output1, "Start\t");}
else{fprintf(output1, "%d\t", count);}
for(k=0;k<num_kins;k++){fprintf(output1, "%.6f\t", hers[k]);}
for(r=0;r<num_regs;r++){fprintf(output1, "%.6f\t", hers[num_kins+r]);}
fprintf(output1, "%f\t", like);
if(count==0){fprintf(output1, "n/a\t%f\t0\n", 0.0001);}
else{fprintf(output1, "%f\t%f\t%d\n", like-like2, 0.0001, nlost);}
fclose(output1);
}

//if(g>0){printf("%f %f\n", hers[num_kins+num_regs], like);}

////////

//now KPY for mkins, regions and gene (probably unnecessary when hers[k]=0)
for(k=0;k<num_kins;k++)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, mkins[k], &ns, PY, &one, &beta, KPY[k], &one);
}
for(r=0;r<num_regs+gflag;r++)
{
token=Xends[r]-Xstarts[r];
XTPY=malloc(sizeof(double)*token);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, PY, &one, &beta, XTPY, &one);
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, XTPY, &one, &beta, KPY[num_kins+r], &one);
free(XTPY);
}

if(adjust>0)	//K[0]PY should become N/N*K[0]PY - N1/N*K[1]PY - ...
{
for(i=0;i<ns;i++)
{
KPY[0][i]=KPY[0][i]*adjust/(adjust-adjust2);
for(r=0;r<num_regs+gflag;r++){KPY[0][i]-=Xsums[r]/(adjust-adjust2)*KPY[1+r][i];}
}
}

////////

//now PKPY for mkins, regions and gene (again, probably unnecessary when hers[k]=0)
for(k=0;k<num_kins+num_regs+gflag;k++)
{
if(shortcut==0)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, KPY[k], &one, &beta, PKPY[k], &one);
}
if(shortcut==1)	//U H UT KPY
{
alpha=1.0;beta=0.0;
if(num_kins==0)	//then U is identity
{
dgemv_("N", &ns, &ns, &alpha, H, &ns, KPY[k], &one, &beta, PKPY[k], &one);
}
else
{
dgemv_("T", &ns, &ns, &alpha, U, &ns, KPY[k], &one, &beta, UTKPY, &one);
dgemv_("N", &ns, &ns, &alpha, H, &ns, UTKPY, &one, &beta, HUTKPY, &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, HUTKPY, &one, &beta, PKPY[k], &one);
}
}
}

////////

//get trace PK for kins, regions and gene (again, unnecessary when hers[k]=0)
for(k=0;k<num_kins;k++)
{
traces[k]=0;
if(shortcut==0)
{
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++){traces[k]+=P[i+i2*ns]*mkins[k][i2+i*ns];}
}
}
if(shortcut==1)	//if here, k=0 so tr(PK) = tr(UHUTUEUT) = tr(HE)
{
for(i=0;i<ns;i++){traces[k]+=E[i]*H[i+i*ns];}
}
}
for(r=0;r<num_regs+gflag;r++)
{
traces[num_kins+r]=0;
token=Xends[r]-Xstarts[r];
PX=malloc(sizeof(double)*ns*token);

if(shortcut==0)	//tr(PK) = tr(PX XT/N)
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &token, &ns, &alpha, P, &ns, X+Xstarts[r]*ns, &ns, &beta, PX, &ns);
for(i=0;i<ns;i++)
{
for(j=0;j<token;j++)
{traces[num_kins+r]+=PX[i+j*ns]*X[i+(Xstarts[r]+j)*ns]/Xsums[r];}
}
}
if(shortcut==1)	//tr(PK) = tr(HUTX XTU) (can store HUTX in PX)
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &token, &ns, &alpha, H, &ns, UTX+Xstarts[r]*ns, &ns, &beta, PX, &ns);
for(i=0;i<ns;i++)
{
for(j=0;j<token;j++)
{traces[num_kins+r]+=PX[i+j*ns]*UTX[i+(Xstarts[r]+j)*ns]/Xsums[r];}
}
}

free(PX);
}

if(adjust>0)	//tr(K[0]P) should become N/N*tr(K[0]P) - N1/N*tr(K[1]P) - ...
{
traces[0]=traces[0]*adjust/(adjust-adjust2);
for(r=0;r<num_regs+gflag;r++){traces[0]-=Xsums[r]/(adjust-adjust2)*traces[1+r];}
}

////////

//add on priors ???

///////////////////////////

//fill AI and BI - set to 0/1 elements corresponding to lambdas fixed at zero
//gam=YTPY, gam2a=YTPVaPY, gam3ab=YTPVaPVbPY
//get AI = -2nd deriv = n/2/gam (-gam2agam2b/gam +gam3ab)

for(k=0;k<num_kins+num_regs+gflag;k++)
{
//diagonals and B
if(herfixes[k]<2)	//else not allowed to update 
{
gam2a=0;for(i=0;i<ns;i++){gam2a+=PKPY[k][i]*Y[i];}
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k][i];}
AI[k+k*(num_kins+num_regs+gflag)]=.5*nfree/gam*(gam3-pow(gam2a,2)/gam);
BI[k]=.5*nfree/gam*gam2a-.5*traces[k];
}
else
{AI[k+k*(num_kins+num_regs+gflag)]=1.0;BI[k]=0.0;}

//off-diagonals
for(k2=0;k2<k;k2++)
{
if(herfixes[k]<2&&herfixes[k2]<2)	//else not allowed to update
{
gam2b=0;for(i=0;i<ns;i++){gam2b+=PKPY[k2][i]*Y[i];}
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k2][i];}
AI[k+k2*(num_kins+num_regs+gflag)]=.5*nfree/gam*(gam3-gam2a*gam2b/gam);
AI[k2+k*(num_kins+num_regs+gflag)]=AI[k+k2*(num_kins+num_regs+gflag)];
}
else
{AI[k+k2*(num_kins+num_regs+gflag)]=0;AI[k2+k*(num_kins+num_regs+gflag)]=0;}
}}	//end of k loop

//save AI but setting zero her terms to zero
for(k=0;k<num_kins+num_regs+gflag;k++)
{
for(k2=0;k2<num_kins+num_regs+gflag;k2++)
{
AIsave[k+k2*(num_kins+num_regs+gflag)]=AI[k2+k*(num_kins+num_regs+gflag)];
if(herfixes[k]>=2){AIsave[k+k2*(num_kins+num_regs+gflag)]=0;}
}
}

//update is invAI BI
(void)eigen_invert(AI, num_kins+num_regs+gflag, AI2, AI3, 0);

token=num_kins+num_regs+gflag;
alpha=1.0;beta=0.0;
dgemv_("N", &token, &token, &alpha, AI, &token, BI, &one, &beta, lamdiffs, &one);

//make sure differences not too large ???
sum=0;for(k=0;k<num_kins+num_regs+gflag;k++){sum+=lambdas[k];}
for(k=0;k<num_kins+num_regs+gflag;k++)
{
if(herfixes[k]<2)
{
if(lamdiffs[k]>.1*pow(1+sum,2)){lamdiffs[k]=.1*pow(1+sum,2);}
if(lamdiffs[k]<-.1*pow(1+sum,2)){lamdiffs[k]=-.1*pow(1+sum,2);}
if(lamdiffs[k]<0.0001-lambdas[k]){lamdiffs[k]=0.0001-lambdas[k];}
lambdas[k]+=lamdiffs[k];
}
}

sum=0;for(k=0;k<num_kins+num_regs+gflag;k++){sum+=lambdas[k];}
for(k=0;k<num_kins+num_regs+gflag;k++)
{
hers[k]=lambdas[k]/(1+sum);

if(lambdas[k]<=0.0001){herfixes[k]++;}
else{herfixes[k]=0;}

if(herfixes[k]==2)	//lambda has been at 0.0001 for two consecutive iterations
{
lambdas[k]=0;hers[k]=0;nlost++;
if(k>=num_kins)	//blank out part of X and for shortcut = 1 corresponding part of UTX
{
for(j=Xstarts[k-num_kins];j<Xends[k-num_kins];j++)
{
for(i=0;i<ns;i++)
{
X[i+j*ns]=0;
if(shortcut==1){UTX[i+j*ns]=0;}
}}
}
}

//get new deltas
if(adjust==0){deltas[k]=lambdas[k];}
if(adjust>0&&k==0){deltas[0]=adjust/(adjust-adjust2)*lambdas[0];}
if(adjust>0&&k>0){deltas[k]=lambdas[k]-Xsums[k-1]/(adjust-adjust2)*lambdas[0];}
}

count++;
}	//end of while loop

if(g==-1)	//screen print
{
printf("Final\t");
for(k=0;k<num_kins;k++){printf("%.6f\t", hers[k]);}
for(r=0;r<num_regs;r++){printf("%.6f\t", hers[num_kins+r]);}
printf("%f\t%f\t%f\t%d\n", like, like-like2, 0.0001, nlost);
}

