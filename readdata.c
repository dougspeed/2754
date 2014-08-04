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


double get_stats_qc(double *centres, double *mults, double *data, int num_samples_use, int length, float minmaf, float maxmaf, float minvar, float minobs, float missingvalue, float power, int ignoremissing, char **prednames, int type, float *weights)
{
//type=1 means standardize, type=2 means standardize with weights
int i, j, indcount;
int misscount, misscount2, misscount3;
double sum, sumsq, mean, var, wsum;
float maf, mindata;

mindata=0;
misscount=0;misscount2=0;misscount3=0;
wsum=0;
for(j=0;j<length;j++)
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

//if requested set missing to mean
if(ignoremissing==1&&mults[j]!=-1&&indcount<num_samples_use)
{
for(i=0;i<num_samples_use;i++)
{
if(data[i+j*num_samples_use]==missingvalue){data[i+j*num_samples_use]=centres[j];}
}
//and update stats
var=sumsq/num_samples_use-pow(mean,2)*indcount/num_samples_use;
}

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

if(type==2)	//standardise and weight
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

