/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//??? want to add in check that fam ids are unique and that snps names are unique?7
//also for kin ids
//improve count els when empty last rows etc


int countrows(char *filename)
{
int count;
char readchar1, readchar2;
FILE *source;

if((source=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

readchar1=0;readchar2=0;
count=0;

while(fscanf(source,"%c",&readchar1)==1)
{
readchar2=readchar1;
if(readchar1==10){count++;}
}
if(readchar2!=10){count++;}	//no end line at last line

fclose(source);

return(count);
}	//end of countrows

/////////

int countels(char *filename)
{
int count;
char readchar1, readchar2;
FILE *source;

if((source=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

readchar1=0;readchar2=0;
count=1;

while(fscanf(source,"%c",&readchar1)==1)
{
if((readchar1==9||readchar1==10||readchar1==32)&&(readchar2!=9&&readchar2!=10&&readchar2!=32)){count++;}
readchar2=readchar1;
}
if(readchar2==9||readchar2==10||readchar2==32){count--;}	//ended in a "space"

fclose(source);

return(count);
}	//end of countels

/////////

int countcols(char *filename)
{
int count;
char readchar1, readchar2;
FILE *source;

if((source=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

readchar1=0;readchar2=0;
count=1;

while(readchar1!=10)
{
readchar2=readchar1;
(void)fscanf(source,"%c", &readchar1);
if((readchar1==9||readchar1==32)&(readchar2!=9&&readchar2!=32)){count++;}
}
if(readchar2==9||readchar2==32){count--;}	//was space after last element

fclose(source);

return(count);
}	//end of countcols

/////////

int read_strings(char *strfile, char **str, int length, int col, int head)
{
int j, k;
char readchar, readstring[100];

FILE *input;

if(countcols(strfile)<col)
{printf("Error reading %s; has %d columns but should have at least %d\n\n", strfile, countels(strfile), col);exit(1);}
if(countrows(strfile)!=length+head)
{printf("Error reading %s; has %d rows but should have %d\n\n", strfile, countrows(strfile), length+head);exit(1);}

if((input=fopen(strfile,"r"))==NULL)
{printf("Error opening %s\n\n", strfile);exit(1);}

if(head==1)
{
readchar=0;while(readchar!=10){(void)fscanf(input, "%c", &readchar);}
}

for(j=0;j<length;j++)
{
for(k=0;k<col-1;k++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s\n\n", k+1, strfile);exit(1);}
}

if(fscanf(input, "%s%c", str[j], &readchar)!=2)
{printf("Error reading element %d of %s\n\n", col, strfile);exit(1);}

//skip to end of row
while(readchar!=10){(void)fscanf(input, "%c", &readchar);}
}

fclose(input);

return(0);
}	//end of read_strings

/////////

int read_ids_merge(char *idfile, char **ids, int length)
{
int i;
char readchar, readstring[100], readstring2[100];

FILE *input;

if(countcols(idfile)<2)
{printf("Error reading %s; has %d columns but should have at least 2\n\n", idfile, countels(idfile));exit(1);}
if(countrows(idfile)!=length)
{printf("Error reading %s; has %d rows but should have %d\n\n", idfile, countrows(idfile), length);exit(1);}

if((input=fopen(idfile,"r"))==NULL)
{printf("Error opening %s\n\n", idfile);exit(1);}

for(i=0;i<length;i++)
{
if(fscanf(input, "%s %s%c", readstring, readstring2, &readchar)!=3)
{printf("Error reading ids from Row %d of %s\n\n", i+1, idfile);exit(1);}
sprintf(ids[i], "%s___%s", readstring, readstring2); 

while(readchar!=10){(void)fscanf(input, "%c", &readchar);}
}

fclose(input);

return(0);
}	//end of read_ids_merge

/////////

int read_keepfile(char *keepfile, int *keepindex, int num_max, int ordered)
{
int j, num_use;
int readint;

FILE *input;

num_use=countels(keepfile);

if((input=fopen(keepfile,"r"))==NULL)
{printf("Error opening file %s\n\n",keepfile);exit(1);}

for(j=0;j<num_use;j++)
{
if(fscanf(input, "%d ", &readint)!=1)
{printf("Error reading element %d of %s\n\n", j+1, keepfile);exit(1);}
if(readint<=0||readint>num_max)
{printf("Error reading %s; element %d has value %d (should be integer between 1 and %d)\n\n", keepfile, j+1, readint, num_max);exit(1);} 
keepindex[j]=readint-1;
}

fclose(input);

if(ordered==1)	//check are ordered
{
for(j=1;j<num_use;j++)
{
if(keepindex[j]<=keepindex[j-1])
{printf("Error reading %s; element %d has value (%d) less than that of the previous element (%d)\n\n", keepfile, j+1, keepindex[j]+1, keepindex[j-1]+1);exit(1);}
}
return(num_use);
}

if(ordered==2)	//make them ordered
{
qsort(keepindex, num_use, sizeof(int), compare_int);
}

return(num_use);
}	//end of read_keepfile

///////////////////////////

int read_famfile (char *famfile,  char **ids1, char **ids2, int num_samples_use, int *keepsamps)
{
int i, count;
char readchar, **ids1b, **ids2b;

FILE *input;


//will read all ids in, then keep the ones we want

count=countrows(famfile);
ids1b=malloc(sizeof(char*)*count);
ids2b=malloc(sizeof(char*)*count);
for(i=0;i<count;i++)
{ids1b[i]=malloc(sizeof(char)*500);ids2b[i]=malloc(sizeof(char)*500);}

if((input=fopen(famfile,"r"))==NULL)
{printf("Error opening %s\n",famfile);exit(1);}

for(i=0;i<count;i++)
{
//read first two elements 
if(fscanf(input, "%s %s%c", ids1b[i], ids2b[i], &readchar)!=3)
{printf("Error reading ids from Row %d of %s\n", i+1, famfile);exit(1);}

//skip to end of line
while(readchar!=10&&i+1!=count){(void)fscanf(input, "%c", &readchar);}
}

fclose(input);

if(keepsamps!=NULL)
{
for(i=0;i<num_samples_use;i++)
{strcpy(ids1[i],ids1b[keepsamps[i]]);strcpy(ids2[i],ids2b[keepsamps[i]]);
}
}
else
{
for(i=0;i<num_samples_use;i++)
{strcpy(ids1[i],ids1b[i]);strcpy(ids2[i],ids2b[i]);}
}

for(i=0;i<count;i++){free(ids1b[i]);free(ids2b[i]);}free(ids1b);free(ids2b);

return(0);
}		//end of read_famfile

/////////

int read_mapfile(char *mapfile, int *chr, char **prednames, double *bp, char *al1, char *al2, int data_length, int *keeppreds_use, int type)
{
//type 0, map; type 1, bim; type 2, map quiet; type 3, bim quiet
int j, count, misscount, misscount2;
char readchar, readchr[100], readbp[100], readstring[100];

int *chrb;
char **prednamesb;
double *bpb;
char *al1b, *al2b;

FILE *input;


//will read in all preds, then extract the ones we want

count=countrows(mapfile);

chrb=malloc(sizeof(int)*count);
prednamesb=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){prednamesb[j]=malloc(sizeof(char)*100);}
bpb=malloc(sizeof(double)*count);
al1b=malloc(sizeof(char)*count);
al2b=malloc(sizeof(char)*count);

if((input=fopen(mapfile,"r"))==NULL)
{printf("Error opening %s\n\n",mapfile);exit(1);}

misscount=0;misscount2=0;
for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s ", readchr, prednamesb[j], readstring)!=3)
{printf("Error reading first three elements of Row %d of %s\n\n", j+1, mapfile);exit(1);}

if(type==0||type==2)	//read bp
{
if(fscanf(input, "%s%c", readbp, &readchar)!=2)
{printf("Error reading basepair of Row %d of %s\n\n", j+1, mapfile);exit(1);}
}

if(type==1||type==3)	//read bp and A1, A2
{
if(fscanf(input, "%s %c %c%c", readbp, al1b+j, al2b+j, &readchar)!=4)
{printf("Error reading basepair and alleles of Row %d of %s\n\n", j+1, mapfile);exit(1);}
}

while(readchar!=10&&j+1!=count)	//skip to end of line
{(void)fscanf(input, "%c", &readchar);}

//save chr
chrb[j]=-1;
if(strcmp(readchr,"X")==0){chrb[j]=23;}
if(strcmp(readchr,"Y")==0){chrb[j]=24;}
if(strcmp(readchr,"XY")==0){chrb[j]=25;}
if(strcmp(readchr,"MT")==0){chrb[j]=26;}
if(strcmp(readchr,"0")==0){chrb[j]=0;}

if(chrb[j]==-1)	//not found so far
{
chrb[j]=atoi(readchr);
if(chrb[j]==0)	//so was not found
{printf("Error reading %s; chromosome provided (%s) on Row %d is not recognised.\nValue must either be positive integer, or X (23), Y (24), XY (25), MT (26) or 0\n\n", mapfile, readchr, j+1);exit(1);}
if(chrb[j]>26&&type<2)	//so unexpectedly large
{
if(misscount<5){printf("Warning reading %s; chromosome provided (%s) on Row %d is greater than 26, the largest expected for humans\n", mapfile, readchr, j+1);}
misscount++;
}
}

//and check it's not smaller than previous one
if(j>0)
{
if(chrb[j]<chrb[j-1])
{printf("Error reading %s, chromosome provided (%s) on Row %d is less than that for previous row (%d)\n", mapfile, readchr, j+1, chrb[j-1]);exit(1);}
}

//save bp and check it's consistent with previous one
bpb[j]=atof(readbp);

if(j>0)
{
if(chrb[j]==chrb[j-1]&&bpb[j]<bpb[j-1])
{printf("Error reading %s; basepair provided (%s) on Row %d is less than that for previous row (%f)\n", mapfile, readbp, j+1, bpb[j-1]);exit(1);}
if(chrb[j]==chrb[j-1]&&bpb[j]==bpb[j-1]&&type<2)
{printf("Warning reading %s; basepair provided (%s) on Row %d is the same as that for previous row\n", mapfile, readbp, j+1);
misscount2++;}
}
}	//end of j loop

fclose(input);

if(misscount>=5){printf("In total %d predictors had chromosome values greater than 26\n", misscount);}
if(misscount>0||misscount2>0){printf("\n");}


for(j=0;j<data_length;j++)
{
chr[j]=chrb[keeppreds_use[j]];strcpy(prednames[j],prednamesb[keeppreds_use[j]]);bp[j]=bpb[keeppreds_use[j]];
if(type!=0){al1[j]=al1b[keeppreds_use[j]];al2[j]=al2b[keeppreds_use[j]];}
}

free(chrb);free(bpb);free(al1b);free(al2b);
for(j=0;j<count;j++){free(prednamesb[j]);};free(prednamesb);

return(0);
}	//end of read_mapfile

/////////

int read_weightfile(char *weightsfile, float *weights, int data_length, int *keeppreds_use, char **prednames, char *mapfile)
{
int j, count, count2;
char readchar, readstring[100];

float *weightsb;

FILE *input;


//will read in all weights, then extract the ones we want

count=countrows(weightsfile);
count2=countcols(weightsfile);

weightsb=malloc(sizeof(float)*count);

if((input=fopen(weightsfile,"r"))==NULL)
{printf("Error opening %s\n",weightsfile);exit(1);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%f%c", weightsb+j, &readchar)!=2)
{printf("Error reading Row %d of %s\n\n", j+1, weightsfile);exit(1);}

if(count2==5)	//read four more elements
{
if(fscanf(input, "%s %s %s %s%c", readstring, readstring, readstring, readstring, &readchar)!=5)
{printf("Error reading last 4 elements of Row %d of %s\n\n", j+1, weightsfile);exit(1);}
if(strcmp(readstring,prednames[j])!=0)
{printf("Error, %s does not align with %s (different SNPs (%s and %s) on Row %d\n\n", weightsfile, mapfile, readstring, prednames[j], j+1);exit(1);}
}
while(readchar!=10){(void)fscanf(input, "%c", &readchar);}
}
fclose(input);


//fill up weights and check none are -1
for(j=0;j<data_length;j++)
{
weights[j]=weightsb[keeppreds_use[j]];
if(weights[j]==-1)
{printf("Error reading %s; the weight for Predictor %s (Row %d) is -1, indicating this predictor was not included when calculating weights, so can not be used in subsequent applications", weightsfile, prednames[j], keeppreds_use[j]+1);}
}

free(weightsb);

return(0);
}	//end of read weightfile

/////////

int read_respfile(char *respfile, double *resp, char **ids1, char **ids2, int num_samples_use, int num_resps_use, int *keepresps, int num_resps, float missingvalue, int **respindex)
{
int i, m, count, found, indcount;
int *keepindex, n0, n1, n2, nm, nr;
double sum, sumsq, mean, var;

char readstring[100];

char **ids1b, **ids2b;
double *resptemp, *respb;

FILE *input;


//will read all ids in, then keep the ones we want

count=countrows(respfile);
ids1b=malloc(sizeof(char*)*count);
ids2b=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){ids1b[i]=malloc(sizeof(char)*100);}
for(i=0;i<count;i++){ids2b[i]=malloc(sizeof(char)*100);}

resptemp=malloc(sizeof(double)*num_resps);
respb=malloc(sizeof(double)*count*num_resps_use);

if((input=fopen(respfile,"r"))==NULL)
{printf("Error opening file %s\n\n",respfile);exit(1);}
printf("Reading phenotypes from %s\n", respfile);

for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", ids1b[i], ids2b[i])!=2)
{printf("Error reading the ids on row %d of %s\n\n", i+1, respfile);exit(1);}

for(m=0;m<num_resps;m++)
{
if(fscanf(input, "%[0-9e.+-] ", readstring)==1)
{resptemp[m]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Phenotype %d for Individual %s %s in %s\n\n", m+1, ids1b[i], ids2b[i], respfile);exit(1);}
if(strcmp(readstring,"NA")==0)
{resptemp[m]=missingvalue;}
else	//it is not NA
{printf("Error reading %s; unknown value (%s) for Phenotype %d, Individual %s %s\n\n", respfile, readstring, m+1, ids1b[i], ids2b[i]);exit(1);}
}
}

//load up resp
for(m=0;m<num_resps_use;m++){respb[i+m*count]=resptemp[keepresps[m]];}
}	//end of i loop

//see which of the ids we can find responses for
keepindex=malloc(sizeof(int)*num_samples_use);
found=find_ids(keepindex, ids1, ids2, num_samples_use, ids1b, ids2b, count);
if(found>0){printf("Error reading %s; does not contain responses for Individual %s %s\n\n", respfile, ids1[found-1], ids2[found-1]);exit(1);}
if(found<0){printf("Error reading %s; features Individual %s %s more than once\n\n", respfile, ids1[-found-1], ids2[-found-1]);exit(1);}

//fill up responses, noting if values are all 0, 1 or 2 (then 0 is missing)
for(m=0;m<num_resps_use;m++)
{
n0=0;n1=0;n2=0;nm=0;
for(i=0;i<num_samples_use;i++)
{
resp[i+m*num_samples_use]=respb[keepindex[i]+m*count];
if(resp[i+m*num_samples_use]==0){n0++;}
if(resp[i+m*num_samples_use]==1){n1++;}
if(resp[i+m*num_samples_use]==2){n2++;}
if(resp[i+m*num_samples_use]==missingvalue){nm++;}
}

nr=num_samples_use-n0-n1-n2-nm;
if(nr==0&&nm==0&&n2>0)
{
printf("Response %d in %s is binary (control:1, case:2, missing:0) - %d controls, %d cases, %d missing\n", keepresps[m]+1, respfile, n1, n2, nm+n0);
for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]==0){resp[i+m*num_samples_use]=missingvalue;}
}
nm+=n0;n0=0;
}
if(nr==0&&n2==0)
{
printf("Response %d in %s is binary (control:0, case:1, missing:%.2f) - %d controls, %d cases, %d missing\n", keepresps[m]+1, respfile, missingvalue, n0, n1, nm);
}
}

//check counts and variances
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue)
{sum+=resp[i+m*num_samples_use];sumsq+=pow(resp[i+m*num_samples_use],2);indcount++;}
}

if(indcount<3)
{printf("Error reading %s; Phenotype %d has only %d non-missing values - come on, you're better than that\n\n", respfile, keepresps[m]+1, indcount);exit(1);}
 
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0)
{printf("Error reading %s; Phenotype %d has variance 0\n\n", respfile, keepresps[m]+1);exit(1);}
}

free(resptemp);free(respb);free(keepindex);
for(i=0;i<count;i++){free(ids1b[i]);free(ids2b[i]);}free(ids1b);free(ids2b);


//record which samples present for each phenotype
for(m=0;m<num_resps_use;m++)
{
count=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue){respindex[m][1+count]=i;count++;}
}
respindex[m][0]=count;
printf("%d out of %d samples have values for Phenotype %d\n", count, num_samples_use, keepresps[m]+1);
}
printf("\n");

return(0);
}	//end of read_respfile

/////////

int read_covarfile(char *covarfile, double *covar, char **ids1, char **ids2, int num_samples_use, int num_covars, float missingvalue)
{
int i, j, count, found, indcount;
int *keepindex;
double sum, sumsq, mean, var;

char readstring[100];

char **ids1b, **ids2b;
double *covarb;

FILE *input;


//will read all ids in, then keep the ones we want

count=countrows(covarfile);
ids1b=malloc(sizeof(char*)*count);
ids2b=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){ids1b[i]=malloc(sizeof(char)*100);}
for(i=0;i<count;i++){ids2b[i]=malloc(sizeof(char)*100);}

covarb=malloc(sizeof(double)*count*num_covars);

if((input=fopen(covarfile,"r"))==NULL)
{printf("Error opening file %s\n\n",covarfile);exit(1);}
printf("Reading covariates from %s\n", covarfile);

for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", ids1b[i], ids2b[i])!=2)
{printf("Error reading the ids on row %d of %s\n\n", i+1, covarfile);exit(1);}

for(j=0;j<num_covars;j++)
{
if(fscanf(input, "%[0-9e.+-] ", readstring)==1)
{covarb[i+j*count]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Covariate %d for Individual %s %s in %s\n\n", j+1, ids1b[i], ids2b[i], covarfile);exit(1);}
if(strcmp(readstring,"NA")==0)
{covarb[i+j*count]=missingvalue;}
else	//it is not NA
{printf("Error reading %s; unknown value (%s) for Covariate %d, Individual %s %s\n\n", covarfile, readstring, j+1, ids1b[i], ids2b[i]);exit(1);}
}
}
}	//end of i loop


//see which of the ids we can find covariates for
keepindex=malloc(sizeof(int)*num_samples_use);
found=find_ids(keepindex, ids1, ids2, num_samples_use, ids1b, ids2b, count);
if(found>0){printf("Error reading %s; does not contain covariates for Individual %s %s\n\n", covarfile, ids1[found-1], ids2[found-1]);exit(1);}
if(found<0){printf("Error reading %s; features Individual %s %s more than once\n\n", covarfile, ids1[-found-1], ids2[-found-1]);exit(1);}

//fill up covariates
for(j=0;j<num_covars;j++)
{
for(i=0;i<num_samples_use;i++){covar[i+j*num_samples_use]=covarb[keepindex[i]+j*count];}
}

//centre and set missing to mean
for(j=0;j<num_covars;j++)
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;
for(i=0;i<num_samples_use;i++)
{
if(covar[i+j*num_samples_use]!=missingvalue)
{sum+=covar[i+j*num_samples_use];sumsq+=pow(covar[i+j*num_samples_use],2);indcount++;}
}

if(indcount<3)
{printf("Error reading %s; Covariate %d has only %d non-missing values - come on, you're better than that\n\n", covarfile, j+1, indcount);exit(1);}
 mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0)
{printf("Error reading %s; Covariate %d has variance 0\n\n", covarfile, j+1);exit(1);}

for(i=0;i<num_samples_use;i++)
{
if(covar[i+j*num_samples_use]!=missingvalue){covar[i+j*num_samples_use]-=mean;}
else{covar[i+j*num_samples_use]=0;}
}
}

free(covarb);free(keepindex);
for(i=0;i<count;i++){free(ids1b[i]);free(ids2b[i]);}free(ids1b);free(ids2b);

return(0);
}	//end of read_covarfile

///////////////////////////

int read_bed_fly(char *bedfile, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, float missingvalue)
{
int i, j, count, rlength, value, value2;
unsigned char check[3], *rowbits;
double conv[4], *datatemp;

FILE *input;


rlength=(int)((num_samples-1)/4)+1;
rowbits=malloc(sizeof(unsigned char)*rlength);

conv[0]=2.0;conv[1]=missingvalue;conv[2]=1.0;conv[3]=0.0;
datatemp=malloc(sizeof(double)*num_samples);

if((input=fopen(bedfile,"rb"))==NULL)
{printf("Error opening file %s\n\n",bedfile);exit(1);}

//check the length and that the first three digits ok
fseek(input, 0, SEEK_END);
if(ftell(input)!=(3+rlength*num_preds)*sizeof(char))
{printf("Error reading %s; does not contain %d x %d SNPs\n\n", bedfile, num_samples, num_preds);exit(1);}
fseek(input, 0, SEEK_SET);
if(fread(check, sizeof(unsigned char), 3, input)!=3)
if(check[0]!=108||check[1]!=27)
{printf("Error reading %s; does not appear to be a binary PLINK file\n\n", bedfile);exit(1);}
if(check[2]!=1)
{printf("Error readinfg %s; can only read in SNP-major mode\n\n", bedfile);exit(1);}

for(j=start;j<end;j++)
{
//if((j-start+1)%10000==0){printf("Reading Predictor %d out of %d from %s\n", j-start+1, end-start, bedfile);}

if(fseek(input, 3+rlength*keeppreds[j], SEEK_SET)!=0)
{printf("Error in %s; unable to find Predictor %d\n\n", bedfile, j+1);exit(1);}
if(fread(rowbits, sizeof(unsigned char), rlength, input)!=rlength)
{printf("Error in %s; unable to read for Predictor %d\n\n", bedfile, j+1);exit(1);}

//convert to basepairs
i=0;
for(count=0;count<rlength;count++)
{
value=(int)rowbits[count];
value2=value%4;
datatemp[i]=conv[value2];i++;
if(i==num_samples){break;}
value2=(value>>2)%4;
datatemp[i]=conv[value2];i++;
if(i==num_samples){break;}
value2=(value>>4)%4;
datatemp[i]=conv[value2];i++;
if(i==num_samples){break;}
value2=(value>>6)%4;
datatemp[i]=conv[value2];i++;
if(i==num_samples){break;}
}

//load up data
for(i=0;i<num_samples_use;i++){data[i+(j-start)*num_samples_use]=datatemp[keepsamps[i]];}
}
fclose(input);

free(rowbits);free(datatemp);

return(0);
}	//end of read_bed_fly

/////////

int read_chiamo_fly(char *chiamofile, FILE *input, int current, int chiamoheaders, int chiamoprobs, double *data, char *al1, char *al2, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, float missingvalue, float threshold)
{
int i, k, num_found, misscount;
float prob0, prob1, prob2, prob3, proball;
double *datatemp;
char readchar, readstring[100];


//will add commands to write to other file types


if(current>keeppreds[start])
{printf("Doug Error reading chiamo - %d past %d\n\n", current, keeppreds[start]);exit(1);}

datatemp=malloc(sizeof(double)*num_samples);

num_found=0;
while(1)
{
while(current<keeppreds[start+num_found])	//skip rows
{
if((current+1)%10000==0){printf("Skipping Row %d of %d of %s\n", current+1, num_preds, chiamofile);}

readchar=0;
while(readchar!=10)
{
if(fscanf(input, "%c", &readchar)!=1)
{printf("Error reading %s; run out of elements!\n\n", chiamofile);exit(1);}
}
current++;
}

//if((current+1)%10000==0){printf("Reading Row %d of %s; %d out of %d predictors kept\n", current+1, chiamofile, num_found, end-start);}

//skip chiamoheaders-2 elements
for(k=0;k<chiamoheaders-2;k++)
{
if(fscanf(input,"%s ", readstring)!=1)
{printf("Error reading first %d headers of Row %d of %s\n\n", current+1, chiamoheaders-2, chiamofile);exit(1);}
}

//now get alleles
if(fscanf(input, "%c %c ", al1+current, al2+current)!=2)
{printf("Error reading the two alleles (cols %d and %d) of %s\n\n", chiamoheaders-1, chiamoheaders, chiamofile);exit(1);}

//now start reading probabilities, and standardise while here
for(i=0;i<num_samples;i++)
{
if(chiamoprobs==2)
{
if(fscanf(input, "%f %f ", &prob0, &prob1)!=2)
{printf("Error in %s reading the two state probabilities for Predictor %d, Individual %d (cols %d and %d)\n\n", chiamofile, current+1, i+1, chiamoheaders+i*2+1, chiamoheaders+i*2+2);exit(1);}
prob2=1-prob0-prob1;
prob3=0;
if(prob0+prob1>1){proball=prob0+prob1;}
else{proball=1;}
}
if(chiamoprobs==3)
{
if(fscanf(input, "%f %f %f ", &prob0, &prob1, &prob2)!=3)
{printf("Error in %s reading the three state probabilities for Predictor %d, Individual %d (cols %d to %d)\n\n", chiamofile, current+1, i+1, chiamoheaders+i*3+1, chiamoheaders+(i+1)*3);exit(1);}
prob3=1-prob0-prob1-prob2;
if(prob0+prob1+prob2>1){proball=prob0+prob1+prob2;}
else{proball=1;}
prob0=prob0/(proball-prob3);
prob1=prob1/(proball-prob3);
prob2=prob2/(proball-prob3);
}
if(chiamoprobs==4)
{
if(fscanf(input, "%f %f %f %f ", &prob0, &prob1, &prob2, &prob3)!=4)
{printf("Error in %s reading the four state probabilities for Predictor %d, Individual %d (cols %d to %d)\n\n", chiamofile, current+1, i+1, chiamoheaders+i*4+1, chiamoheaders+(i+1)*4);exit(1);}
proball=prob0+prob1+prob2+prob3;
prob0=prob0/(proball-prob3);
prob1=prob1/(proball-prob3);
prob2=prob2/(proball-prob3);
}

if(proball>1.01)
{printf("Error reading %s; probabilities for Predictor %d, Individual %d (cols %d to %d) sum to %f (greater than 1)\n\n", chiamofile, current+1, i+1, chiamoheaders+i*chiamoprobs+1, chiamoheaders+(i+1)*chiamoprobs, proball);exit(1);}

if(proball-prob3<0.98&&threshold==-1)
{
if(misscount<10)
{printf("Error reading %s; state probabilities for Predictor %d, Individual %d (cols %d to %d) sum to only %f (less than 0.98) so will treat as missing. This is uncommon using imputed data.\n", chiamofile, current+1, i+1, chiamoheaders+i*chiamoprobs+1, chiamoheaders+(i+1)*chiamoprobs, proball);}
misscount++;
}

datatemp[i]=prob1+2*prob0;
if(proball-prob3<0.98){datatemp[i]=missingvalue;}

if(threshold!=-1)	//want hard calls
{
if(prob0>=threshold){datatemp[i]=0;}
if(prob1>=threshold){datatemp[i]=1;}
if(prob2>=threshold){datatemp[i]=2;}
}
}

//load up data
for(i=0;i<num_samples_use;i++){data[i+(num_found)*num_samples_use]=datatemp[keepsamps[i]];}

num_found++;current++;
if(num_found==end-start){break;}
}

free(datatemp);

return(current);
}	//end of read_chiamo_fly

/////////

int read_sp_fly(char *spfile, FILE *input, int current, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, float missingvalue)
{
int i, num_found;
double *datatemp;
char readchar, readstring[100];


//will add commands to write to other file types


if(current>keeppreds[start])
{printf("Doug Error reading sp - %d past %d\n\n", current, keeppreds[start]);exit(1);}

datatemp=malloc(sizeof(double)*num_samples);

num_found=0;
while(1)
{
while(current<keeppreds[start+num_found])	//skip rows
{
if((current+1)%10000==0){printf("Skipping Row %d of %d of %s\n", current+1, num_preds, spfile);}

readchar=0;
while(readchar!=10)
{
if(fscanf(input, "%c", &readchar)!=1)
{printf("Error reading %s; run out of elements!\n\n", spfile);exit(1);}
}
current++;
}


for(i=0;i<num_samples;i++)
{
if(fscanf(input, "%[0-9e.+-] ", readstring)==1)	
{datatemp[i]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading value for Predictor %d, Individual %d in %s\n\n", current+1, i+1, spfile);}
if(strcmp(readstring,"NA")==0)
{datatemp[i]=missingvalue;}
else	//it is not NA
{printf("Error reading %s; unknown value (%s) for Predictor %d, Individual %d\n", spfile, readstring, current+1, i+1);exit(1);}
}
}

//load up data
for(i=0;i<num_samples_use;i++){data[i+(num_found)*num_samples_use]=datatemp[keepsamps[i]];}

num_found++;current++;
if(num_found==end-start){break;}
}

free(datatemp);

return(current);
}	//end of read_chiamo_fly

/////////

int read_speed_fly(char *speedfile, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, float missingvalue)
{
int i, j;
double *datatemp;

FILE *input;

datatemp=malloc(sizeof(double)*num_samples);

if((input=fopen(speedfile,"rb"))==NULL)
{printf("Error opening file %s\n\n",speedfile);exit(1);}

//check the length
fseek(input, 0, SEEK_END);
if(ftell(input)!=num_samples*num_preds*sizeof(double))
{printf("Error reading %s; does not contain %d x %d SNPs\n\n", speedfile, num_samples, num_preds);exit(1);}
fseek(input, 0, SEEK_SET);

for(j=start;j<end;j++)
{
//if((j-start+1)%10000==0){printf("Reading Predictor %d out of %d from %s\n", j-start+1, end-start, speedfile);}

if(fseek(input, num_samples*keeppreds[j]*sizeof(double), SEEK_SET)!=0)
{printf("Error in %s; unable to find Predictor %d\n\n", speedfile, j+1);exit(1);}
if(fread(datatemp, sizeof(double), num_samples, input)!=num_samples)
{printf("Error in %s; unable to read for Predictor %d\n\n", speedfile, j+1);exit(1);}

//load up data
for(i=0;i<num_samples_use;i++){data[i+(j-start)*num_samples_use]=datatemp[keepsamps[i]];}
}
fclose(input);

free(datatemp);

return(0);
}	//end of read_speed_fly

///////////////////////////

int read_genefile(char *genefile, char **genenames, int *gchr, double *gbp1, double *gbp2, int num_genes)
{
int j;
int count, misscount;

char readchar, readchr[100], readbp1[100], readbp2[100];

FILE *input;

//no need to check number of lines, as num_genes obtained from this number

if((input=fopen(genefile,"r"))==NULL)
{printf("Error opening file %s\n\n",genefile);exit(1);}

count=0;misscount=0;
for(j=0;j<num_genes;j++)
{
//read the gene name
if(fscanf(input, "%s ", genenames[j])!=1)
{printf("Error reading the gene name for row %d of %s\n\n", j+1, genefile);exit(1);}

//read the chromosome number
if(fscanf(input, "%s ", readchr)!=1)
{printf("Error reading the chromosome for row %d of %s (gene %s)\n\n", j+1, genefile, genenames[j]);exit(1);}

//read the basepairs
if(fscanf(input, "%s %s%c", readbp1, readbp2, &readchar)!=3)
{printf("Error reading the start and end basepairs for row %d of %s (gene %s)\n\n", j+1, genefile, genenames[j]);exit(1);}

while(readchar!=10&&j+1!=num_genes)	//skip to end of line
{(void)fscanf(input, "%c", &readchar);}

//save chr
gchr[j]=0-1;
if(strcmp(readchr,"X")==0){gchr[j]=23;}
if(strcmp(readchr,"Y")==0){gchr[j]=24;}
if(strcmp(readchr,"XY")==0){gchr[j]=25;}
if(strcmp(readchr,"MT")==0){gchr[j]=26;}
if(strcmp(readchr,"0")==0){gchr[j]=0;}

if(gchr[j]==-1)	//not found so far
{
gchr[j]=atoi(readchr);
if(gchr[j]==0)	//so was not found
{printf("Error reading %s; chromosome provided (%s) on Row %d is not recognised.\nValue must either be positive integer, or X (23), Y (24), XY (25), MT (26) or 0\n\n", genefile, readchr, j+1);exit(1);}
if(gchr[j]>26&&misscount<5)	//so unexpectedly large
{
printf("Warning reading %s; chromosome provided (%s) on Row %d is greater than 26, the largest expected for humans\n", genefile, readchr, j+1);
misscount++;
}
}

//and check it's not smaller than previous one
if(j>0)
{
if(gchr[j]<gchr[j-1])
{printf("Error reading %s; chromosome provided (%s) for Gene %s (Row %d) is less than that for previous gene (%d)\n", genefile, readchr, genenames[j], j+1, gchr[j-1]);exit(1);}
}

//save bps
gbp1[j]=atof(readbp1);
gbp2[j]=atof(readbp2);

if(j>0)	//check bp1 is not smaller than previous one
{
if(gchr[j]==gchr[j-1]&&gbp1[j]<gbp1[j-1])
{printf("Error reading %s; base pair provided (%s) for Gene %s (Row %d) is less than that for previous gene (%f)\n\n", genefile, readbp1, genenames[j], j+1, gbp1[j-1]);exit(1);}
}

//and that bp2 is not smaller than bp1
if(gbp2[j]<gbp1[j])
{printf("Error reading %s; start base pair (%s) provided for Gene %s (Row %d) is greater than end base pair (%s)\n\n", genefile, readbp1, genenames[j], j+1, readbp2);exit(1);}

if(gbp1[j]==gbp2[j]&&count<10)
{printf("Warning reading %s; Gene %s (Row %d) has length zero (both start and end base pairs are %s) - I hope this is correct!\n", genefile, genenames[j], j+1, readbp1);count++;}

}	//end of j loop

fclose(input);

return(0);
}	//end of read_genefile

////////

int read_kinfile(char *kinstem, double *kins, double *divs, int num_samples_use, char **ids1, char **ids2, int type)
//type=0, just need kins; type=1, add; type=2, subtract
{
int i, iwant, i2, count, count2;
int *kindex, *kindexb;

float *datatemp1, *datatemp2;
char **usesamps, **allsamps;

char filename1[500], filename2[500], filename3[500];
FILE *input2, *input3;


//first work out which individuals we want
usesamps=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++){usesamps[i]=malloc(sizeof(char)*200);}
for(i=0;i<num_samples_use;i++){sprintf(usesamps[i], "%s___%s", ids1[i], ids2[i]);}

sprintf(filename1, "%s.grm.id", kinstem);
count=countrows(filename1);
allsamps=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){allsamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(filename1, allsamps, count);

kindex=malloc(sizeof(int)*count);
kindexb=malloc(sizeof(int)*count);
count2=find_strings(allsamps, count, usesamps, num_samples_use, kindex, kindexb);
if(count2!=num_samples_use)
{printf("Error reading %s; only found %d of the %d individuals\n\n", filename1, count2, num_samples_use);exit(1);}

/////////

datatemp1=malloc(sizeof(float)*count);
datatemp2=malloc(sizeof(float)*count);

//open and check size of kinship file
sprintf(filename2, "%s.grm.bin", kinstem);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening file %s\n\n", filename2);exit(1);}
fseek(input2, 0, SEEK_END);
if(ftell(input2)!=count*(count+1)/2*sizeof(float))
{printf("Error reading %s; does not contain kinships for %d samples\n\n", filename2, count);exit(1);}
fseek(input2, 0, SEEK_SET);

if(type==1||type==2)	//and also of counts file
{
sprintf(filename3, "%s.grm.N.bin", kinstem);
if((input3=fopen(filename3,"rb"))==NULL)
{printf("Error opening file %s\n\n", filename3);exit(1);}
fseek(input3, 0, SEEK_END);
if(ftell(input3)!=count*(count+1)/2*sizeof(float))
{printf("Error reading %s; does not contain counts for %d samples\n\n", filename3, count);exit(1);}
fseek(input3, 0, SEEK_SET);
}

for(i=0;i<num_samples_use;i++)
{
iwant=kindex[i];

fseek(input2, iwant*(iwant+1)/2*sizeof(float), SEEK_SET);
if(fread(datatemp1, sizeof(float), iwant+1, input2)!=iwant+1)
{printf("Doug error reading 1\n");exit(1);}
if(type==1||type==2)
{
fseek(input3, iwant*(iwant+1)/2*sizeof(float), SEEK_SET);
if(fread(datatemp2, sizeof(float), iwant+1, input3)!=iwant+1)
{printf("Doug error reading 2\n");exit(1);}
}

if(type==0)	//just need kins
{
for(i2=0;i2<=i;i2++){kins[kindexb[i]+kindexb[i2]*num_samples_use]=datatemp1[kindex[i2]];}
for(i2=0;i2<i;i2++){kins[kindexb[i2]+kindexb[i]*num_samples_use]=datatemp1[kindex[i2]];}
}
if(type==1)	//add kins and counts
{
for(i2=0;i2<=i;i2++)
{kins[kindexb[i]+kindexb[i2]*num_samples_use]+=datatemp1[kindex[i2]]*datatemp2[kindex[i2]];
divs[kindexb[i]+kindexb[i2]*num_samples_use]+=datatemp2[kindex[i2]];}
for(i2=0;i2<i;i2++)
{kins[kindexb[i2]+kindexb[i]*num_samples_use]+=datatemp1[kindex[i2]]*datatemp2[kindex[i2]];
divs[kindexb[i2]+kindexb[i]*num_samples_use]+=datatemp2[kindex[i2]];}
}
if(type==2)	//subtract kins and counts
{
for(i2=0;i2<=i;i2++)
{kins[kindexb[i]+kindexb[i2]*num_samples_use]-=datatemp1[kindex[i2]]*datatemp2[kindex[i2]];
divs[kindexb[i]+kindexb[i2]*num_samples_use]-=datatemp2[kindex[i2]];}
for(i2=0;i2<i;i2++)
{kins[kindexb[i2]+kindexb[i]*num_samples_use]-=datatemp1[kindex[i2]]*datatemp2[kindex[i2]];
divs[kindexb[i2]+kindexb[i]*num_samples_use]-=datatemp2[kindex[i2]];}
}
}

fclose(input2);
if(type==1||type==2){fclose(input3);}

for(i=0;i<num_samples_use;i++){free(usesamps[i]);}free(usesamps);
for(i=0;i<count;i++){free(allsamps[i]);}free(allsamps);
free(kindex);free(kindexb);free(datatemp1);free(datatemp2);

return(0);
}	//end of read_kinfile

/////////

int read_partition(char *partfile, int *pindex, char **prednames, int num_preds_use, int partition)
{
int j, count, count2;

char **usesnps;


count=countrows(partfile);
usesnps=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){usesnps[j]=malloc(sizeof(char)*100);}
read_strings(partfile, usesnps, count, 1, 0);

count2=find_strings(prednames, num_preds_use, usesnps, count, pindex, NULL);
if(count2<count)
{printf("Error, only found %d of the %d predictor ids listed in %s\n\n", count2, count, partfile);exit(1);}
printf("Partition list %d, stored in %s, contains %d predictors\n", partition, partfile, count);

for(j=0;j<count;j++){free(usesnps[j]);}free(usesnps);

return(count);
}	//end of read_partition

/////////

int read_regions(char *regpref, int num_regs, char *mapfile, int **regindex, int *rkeeppreds_use)
{
int j, r, length1, length2, max, found, count;
int *usedpreds, *rev;
char **allsnps, **usesnps;

char filename[500];

length1=countrows(mapfile);
allsnps=malloc(sizeof(char*)*length1);
for(j=0;j<length1;j++){allsnps[j]=malloc(sizeof(char)*100);}
read_strings(mapfile, allsnps, length1, 2, 0);

max=0;
for(r=0;r<num_regs;r++)
{
sprintf(filename, "%s%d", regpref, r+1);
length2=countrows(filename);
if(length2>max){max=length2;}
usesnps=malloc(sizeof(char*)*length2);
for(j=0;j<length2;j++){usesnps[j]=malloc(sizeof(char)*100);}
read_strings(filename, usesnps, length2, 1, 0);

regindex[r]=malloc(sizeof(int)*(1+length2));
regindex[r][0]=find_strings(allsnps, length1, usesnps, length2, regindex[r]+1, NULL);
if(regindex[r][0]<length2)
{printf("Error, only found %d of the %d predictors stored in %s\n\n", regindex[r][0], length2, filename);exit(1);}
printf("Region list %d, stored in %s, contains %d predictors\n", r+1, filename, length2);

for(j=0;j<length2;j++){free(usesnps[j]);}free(usesnps);
}
printf("\n");

//tally which ones used and get reverse mapping
usedpreds=malloc(sizeof(int)*length1);
for(j=0;j<length1;j++){usedpreds[j]=0;}
rev=malloc(sizeof(int)*length1);

for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++){usedpreds[regindex[r][1+j]]++;}
}

found=0;count=0;
for(j=0;j<length1;j++)
{
if(usedpreds[j]>0){rkeeppreds_use[found]=j;rev[j]=found;found++;}
if(usedpreds[j]>1)
{
if(count<10)
{printf("Warning, Predictor %s specified in more than one region list\n", allsnps[j]);}
count++;
}
}

if(count>0)
{printf("\nIn total, %d predictors were specified in more than one region list; when specifying regions using the arguments \"--region-number\" and \"--region-prefix\", these are usually (but not necessarily) non-overlapping.\n\n", count);}

//make regindexes relative to rkeeppreds_use
for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++){regindex[r][1+j]=rev[regindex[r][1+j]];}
}

for(j=0;j<length1;j++){free(allsnps[j]);}free(allsnps);
free(usedpreds);free(rev);

return(found);
}	//end of read_regions

/////////

int read_eigenfile(char *eigenfile, double *U, double *E, int num_samples_use, char **ids1, char **ids2)
{
int i, i2;
float readfloat;
char readid1[100], readid2[100];

//double *Ub, *Eb;
//char **ids1b, **ids2b;

FILE *input;

//need to update this to allow diff orders


if(countcols(eigenfile)!=num_samples_use+3)
{printf("Error, saved eigen-decomposition (%s) should have %d columns (not %d) - are you trying to use a different phenotype file?\n\n", eigenfile, num_samples_use+3, countcols(eigenfile));exit(1);}

/*
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);
ids1b=malloc(sizeof(char*)*num_samples_use);
ids2b=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{ids1b[i]=malloc(sizeof(char)*500);ids2b[i]=malloc(sizeof(char)*500);}
*/

if((input=fopen(eigenfile,"r"))==NULL)
{printf("Error opening %s\n\n", eigenfile);exit(1);}

for(i=0;i<num_samples_use;i++)
{
if(fscanf(input, "%s %s ", readid1, readid2)!=2)
{printf("Error reading ids from Row %d of %s\n\n", i+1, eigenfile);exit(1);}
if(strcmp(ids1[i], readid1)!=0||strcmp(ids1[i], readid1)!=0)
{printf("Error reading %s; ids on Row %d should equal %s %s (not %s %s)\n\n", eigenfile, i+1, ids1[i], ids2[i], readid1, readid2);exit(1);}

if(fscanf(input, "%f ", &readfloat)!=1)
{printf("Error reading Row %d Column %d of %s\n\n", i+1, i2+3, eigenfile);exit(1);}
E[i]=readfloat;

for(i2=0;i2<num_samples_use;i2++)
{
if(fscanf(input, "%f ", &readfloat)!=1)
{printf("Error reading Row %d Column %d of %s\n\n", i+1, i2+4, eigenfile);exit(1);}
U[i+i2*num_samples_use]=readfloat;
}

}
fclose(input);

return(0);
}	//end of read_eigenfile


