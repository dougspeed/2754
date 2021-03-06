/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/



int sort_coding(double *data, int num_samples_use, int data_length, int encoding, float missingvalue)
{
int i, j;

if(encoding==2)	//convert 0/1/2 to 0/2/2
{
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<data_length;j++)
{
if(data[i+j*num_samples_use]!=0&&data[i+j*num_samples_use]!=1&&data[i+j*num_samples_use]!=2&&data[i+j*num_samples_use]!=missingvalue){printf("Error recoding to dominant model. Predictor %d for individual %d takes value %f (should be 0, 1 or 2)\n\n", j+1, i+1, data[i+j*num_samples_use]);exit(1);}
if(data[i+j*num_samples_use]==1){data[i+j*num_samples_use]=2;}
}
}}

if(encoding==3)	//convert 0/1/2 to 0/0/2
{
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<data_length;j++)
{
if(data[i+j*num_samples_use]!=0&&data[i+j*num_samples_use]!=1&&data[i+j*num_samples_use]!=2&&data[i+j*num_samples_use]!=missingvalue){printf("Error recoding to recessive model. Predictor %d for individual %d takes value %f (should be 0, 1 or 2)\n\n", j+1, i+1, data[i+j*num_samples_use]);exit(1);}
if(data[i+j*num_samples_use]==1){data[i+j*num_samples_use]=0;}
}
}}

if(encoding==4)	//convert 0/1/2 to 0/2/0
{
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<data_length;j++)
{
if(data[i+j*num_samples_use]!=0&&data[i+j*num_samples_use]!=1&&data[i+j*num_samples_use]!=2&&data[i+j*num_samples_use]!=missingvalue){printf("Error recoding to heterozygote model. Predictor %d for individual %d takes value %f (should be 0, 1 or 2)\n\n", j+1, i+1, data[i+j*num_samples_use]);exit(1);}
if(data[i+j*num_samples_use]==2){data[i+j*num_samples_use]=0;}
if(data[i+j*num_samples_use]==1){data[i+j*num_samples_use]=2;}
}
}}

return(0);
}	//end of sort_coding

///////////////////////////

double get_stats_qc(double *centres, double *mults, double *data, int num_samples_use, int length, float minmaf, float maxmaf, float minvar, float minobs, float missingvalue, float power, char **prednames, int type, float *weights)
{
//type=0, just get stats, type=1 means standardize, type=2 means standardize with weights, type=3, standardize with weights using orig stats
int i, j, indcount;
int misscount, misscount2, misscount3;
double sum, sumsq, mean, var, wsum;
float maf, mindata;

mindata=0;
misscount=0;misscount2=0;misscount3=0;
wsum=0;
for(j=0;j<length;j++)
{
if(type==0||type==1||type==2)	//will get stats
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;maf=0;
for(i=0;i<num_samples_use;i++)
{
if(data[i+j*num_samples_use]!=missingvalue)
{
sum+=data[i+j*num_samples_use];sumsq+=pow(data[i+j*num_samples_use],2);indcount++;
if(data[i+j*num_samples_use]<mindata){mindata=data[i+j*num_samples_use];}
}
}
if(indcount>0)
{
mean=sum/indcount;
maf=mean/2;if(maf>0.5&&maf<=1){maf=1-maf;}
if(minmaf==0)	//make sure not filtered based on maf
{maf=maxmaf/2;}
}
if(indcount>0){var=sumsq/indcount-pow(mean,2);}

//set centres and mults, but can still fail in next step
centres[j]=mean;
mults[j]=pow(var,power/2);

if(indcount<(int)(minobs*num_samples_use))	//too many missing
{
if(misscount<10)
{printf("Warning, Predictor %s has too many missing - will be ignored\n", prednames[j]);}
misscount++;
mults[j]=-1;
}

if(var==0)	//trivial
{
if(misscount2<10&&mults[j]!=-1)
{printf("Warning, Predictor %s is trivial (takes at most one non-missing value across samples) - will be ignored\n", prednames[j]);
}
if(mults[j]!=-1){misscount2++;}
mults[j]=-1;
}

if(var<minvar||maf<minmaf||maf>maxmaf)	//fail QC
{
if(misscount3<10&&mults[j]!=-1)
{printf("Warning, Predictor %s shows insufficient variation (MAF %f, Variance %f) - will be ignored\n", prednames[j], maf, var);}
if(mults[j]!=-1){misscount3++;}
mults[j]=-1;
}

//will set missing to mean so update stats
if(mults[j]!=-1)
{var=sumsq/num_samples_use-pow(mean,2)*indcount/num_samples_use;
mults[j]=pow(var,power/2);}
}	//end of getting stats

if(type==1)	//standardise
{
if(mults[j]!=-1)
{
for(i=0;i<num_samples_use;i++)
{
if(data[i+j*num_samples_use]!=missingvalue)
{data[i+j*num_samples_use]=(data[i+j*num_samples_use]-centres[j])*mults[j];}
else
{data[i+j*num_samples_use]=0;}
}
}
else
{
for(i=0;i<num_samples_use;i++){data[i+j*num_samples_use]=0;}
}
}

if(type==2||type==3)	//standardise and weight
{
if(mults[j]!=-1&&weights[j]>0)
{
wsum+=weights[j];
for(i=0;i<num_samples_use;i++)
{
if(data[i+j*num_samples_use]!=missingvalue)
{data[i+j*num_samples_use]=(data[i+j*num_samples_use]-centres[j])*mults[j]*pow(weights[j],.5);}
else
{data[i+j*num_samples_use]=0;}
}
}
else
{
for(i=0;i<num_samples_use;i++){data[i+j*num_samples_use]=0;}
}
}
}	//end of j loop

return(wsum);
}	//end of get_stats_qc

///////////////////////////

FILE* open_data_fly(char *chiamofile, char *spfile, int num_samples, int chiamoprobs, int chiamoheaders)
{
FILE *input;

if(strcmp(chiamofile,"blank")!=0)
{
if(countcols(chiamofile)!=num_samples*chiamoprobs+chiamoheaders)
{printf("Error reading %s; should contain %d columns, not %d\n\n", chiamofile, num_samples*chiamoprobs+chiamoheaders, countcols(chiamofile));exit(1);}

if((input=fopen(chiamofile,"r"))==NULL)
{printf("Error opening %s\n",chiamofile);exit(1);}
}

if(strcmp(spfile,"blank")!=0)
{
if(countcols(spfile)!=num_samples)
{printf("Error reading %s; should contain %d columns, not %d\n\n", spfile, num_samples, countcols(spfile));exit(1);}

if((input=fopen(spfile,"r"))==NULL)
{printf("Error opening %s\n",spfile);exit(1);}
}

return(input);
}

///////////////////////////

int read_data_fly(char *bedfile, char *chiamofile, char *spfile, char *speedfile, FILE *datainput, int current, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds_use, int num_samples, int num_preds, float missingvalue, int chiamoheaders, int chiamoprobs, char *al1, char *al2, float threshold)
{
if(strcmp(bedfile,"blank")!=0)
{read_bed_fly(bedfile, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue);}
if(strcmp(chiamofile,"blank")!=0)
{
current=read_chiamo_fly(chiamofile, datainput, current, chiamoheaders, chiamoprobs, data, al1, al2, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue, threshold);
}
if(strcmp(spfile,"blank")!=0)
{
current=read_sp_fly(spfile, datainput, current, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue);
}
if(strcmp(speedfile,"blank")!=0)
{read_speed_fly(speedfile, data, num_samples_use, keepsamps, start, end, keeppreds_use, num_samples, num_preds, missingvalue);}

return(current);
}	//end of read_data_fly

///////////////////////////

double cor_pair(double *data, int j, int k, int num_samples_use)
{
int i;
double sumj, sumk, sumjk, sumjj, sumkk;
double meanj, meank, varj, vark, cor;

sumj=0;sumk=0;sumjk=0;sumjj=0;sumkk=0;
for(i=0;i<num_samples_use;i++)
{
sumj+=data[i+j*num_samples_use];
sumk+=data[i+k*num_samples_use];
sumjk+=data[i+j*num_samples_use]*data[i+k*num_samples_use];
sumjj+=pow(data[i+j*num_samples_use],2);
sumkk+=pow(data[i+j*num_samples_use],2);
}
meanj=sumj/num_samples_use;meank=sumk/num_samples_use;
varj=sumjj-sumj*meanj;vark=sumkk-sumk*meank;
cor=pow(sumjk-sumj*meank,2)/varj/vark;

return(cor);
}

/////////

int prune_regions(int num_regs, int **regindex, double *rdata, int num_samples_use, int rdata_length, float prune, double *rmults, float *rweights)
{
//prune and change weights, but can't change weight of snps appearing in more than one region
int j, j2, k, k2, r, count, count2, count3;
int found, keep2, lose, lose2;
int *usedpreds, *retain;

usedpreds=malloc(sizeof(int)*rdata_length);
for(j=0;j<rdata_length;j++){usedpreds[j]=0;}
for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++){usedpreds[regindex[r][1+j]]++;}
}

if(prune<=1){printf("Pruning SNPs based on a correlation squared of %f\nTo avoid this, set \"--prune\" to a value greater than 1\n", prune);}

count=0;count3=0;
for(r=0;r<num_regs;r++)	//loop through regions
{
count3+=regindex[r][0];
retain=malloc(sizeof(int)*regindex[r][0]);
for(j=0;j<regindex[r][0];j++)
{
retain[j]=1;
j2=regindex[r][1+j];
if(rweights[j2]==0||rmults[j2]==-1){retain[j]=0;}
}

if(prune<=1)
{
printf("Scanning through Region %d\n", r+1);
for(j=0;j<regindex[r][0];j++)
{
if(retain[j]==1)
{
j2=regindex[r][1+j];
for(k=0;k<j;k++)
{
if(retain[k]==1)
{
k2=regindex[r][1+k];

//get correlation - will be no missing
if(cor_pair(rdata, j2, k2, num_samples_use)>prune)
{
//can we keep either of j2 or k2 (ie used only once)?
found=0;
if(usedpreds[j2]==1){found=1;keep2=j2;lose=k;lose2=k2;}
if(usedpreds[k2]==1){found=1;keep2=k2;lose=j;lose2=j2;}

if(found==1){rweights[keep2]+=rweights[lose2];rweights[lose2]=0;retain[lose]=0;usedpreds[lose2]--;}
}
}}	//end of retain[k]=1 and k loop}
}}	//end of retain[j]=1 and j loop
}	//end of prune<=1

//refill region
count2=0;
for(j=0;j<regindex[r][0];j++)
{
if(retain[j]==1){regindex[r][1+count2]=regindex[r][1+j];count2++;}
}
regindex[r][0]=count2;
count+=count2;
free(retain);
}	//end of r loop

if(prune<=1){printf("Number of region predictors reduced from %d to %d\n\n", count3, count);}
free(usedpreds);

return(0);
}	//end of prune_regions

///////////////////////////

int fill_X(double *X, int *Xstarts, int *Xends, int *Xrec, int *Xrev, double *Xsums, int *kindex, int num_samples_use, int num_regs, double *rdata, int **regindex, double *rmults, float *rweights, double *gdata, int glength, double *gmults, float *gweights, float *adjust2)
{
int i,j, j2, r, count;
int ns, wnum;

ns=kindex[0];
*adjust2=0;

wnum=0;*adjust2=0;count=0;
for(r=0;r<num_regs;r++)
{
if(Xstarts!=NULL){Xstarts[r]=wnum;Xsums[r]=0;}
for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
if(rmults[j2]!=-1&&rweights[j2]>0)
{
for(i=0;i<ns;i++){X[i+wnum*ns]=rdata[kindex[1+i]+j2*ns];}
if(Xstarts!=NULL)
{Xrec[wnum]=r;Xrev[wnum]=j2;Xsums[r]+=rweights[j2];}
*adjust2+=rweights[j2];
wnum++;
}}
if(Xstarts!=NULL){Xends[r]=wnum;}
}
if(glength>0)
{
if(Xstarts!=NULL){Xstarts[num_regs]=wnum;Xsums[num_regs]=0;}
for(j=0;j<glength;j++)
{
if(gmults[j]!=-1&&gweights[j]>0)
{
for(i=0;i<ns;i++){X[i+wnum*ns]=gdata[kindex[1+i]+j*ns];}
if(Xstarts!=NULL)
{Xrec[wnum]=num_regs;Xrev[wnum]=j;Xsums[num_regs]+=gweights[j];}
*adjust2+=gweights[j];
wnum++;count++;
}}
if(Xstarts!=NULL){Xends[num_regs]=wnum;}
}

return(count);
}	//end of fill_X


