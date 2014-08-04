/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int set_up_blup1(float **allcentres, float **allfactors, float *wsums, int num_kins, char **kinstems, int num_regs, char *regfile, int num_preds_use, int *keeppreds, int *keeppreds_use, char **prednames, char *al1, char *al2, char *famfile, float adjust)
{
int j, j2, j3, k, r, count, count2, flag;
int *kindex, *kindexb, *usedpreds;

int readint;
char readchar, readstring[100];

float *ctemp, *mtemp, *wtemp, *ftemp;
char **ptemp, *a1temp, *a2temp;

FILE *input1, *input2;
char filename1[500];


//read in details, then see which ones we will be using
usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

for(k=0;k<num_kins+num_regs;k++)
{
allcentres[k]=malloc(sizeof(float)*num_preds_use);
allfactors[k]=malloc(sizeof(float)*num_preds_use);
for(j=0;j<num_preds_use;j++){allcentres[k][j]=0;allfactors[k][j]=0;}
}

for(k=0;k<num_kins;k++)	//read in each kinship file
{
sprintf(filename1, "%s.grm.details", kinstems[k]);
printf("Reading details for Kinship %d from %s\n", k+1, filename1);
count=countrows(filename1)-1;

ptemp=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){ptemp[j]=malloc(sizeof(char)*100);}
ctemp=malloc(sizeof(float)*count);mtemp=malloc(sizeof(float)*count);
wtemp=malloc(sizeof(float)*count);
a1temp=malloc(sizeof(char)*count);a2temp=malloc(sizeof(char)*count);

kindex=malloc(sizeof(int)*count);kindexb=malloc(sizeof(int)*count);

//open and skip header line
if((input1=fopen(filename1,"r"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input1,"%c", &readchar);}

flag=0;count2=0;
for(j=0;j<count;j++)
{
if(fscanf(input1, "%s %d %f %f %f %c %c ", ptemp[j], &readint, ctemp+j, mtemp+j, wtemp+j, a1temp+j, a2temp+j)!=7)
{printf("Error reading Row %d of %s\n\n", j+2, filename1);exit(1);}
readint--;

if(a1temp[j]==a2temp[j])
{printf("Error reading %s, both alleles for %s are the same (%c)\n\n", filename1, ptemp[j], a1temp[j]);exit(1);}

if(flag==0&&readint>=0&&readint<num_preds_use)
{
if(strcmp(ptemp[j], prednames[readint])==0)
{kindex[count2]=j;kindexb[count2]=readint;count2++;}
else
{flag=1;}
}
else
{flag=1;}
}
fclose(input1);

if(flag==1)	//could not find predictors using index, so find manually
{
count2=find_strings(ptemp, count, prednames, num_preds_use, kindex, kindexb);
if(count2<count)
{printf("Data are available for %d of the %d predictors\n", count2, count);}
}	//end of flag=1

//now check valid and load up
wsums[k]=0;
for(j=0;j<count;j++)
{
j2=kindex[j];j3=kindexb[j];
if(a1temp[j2]!=al1[j3]&&a1temp[j2]!=al2[j3])
{printf("Error reading %s, alleles for %s (%c and %c) not consistent with those in data (%c and %c)\n\n", filename1, ptemp[j2], a1temp[j2], a2temp[j2], al1[j3], al2[j3]);exit(1);}
if(a2temp[j2]!=al1[j3]&&a2temp[j2]!=al2[j3])
{printf("Error reading %s, alleles for %s (%c and %c) not consistent with those in data (%c and %c)\n\n", filename1, ptemp[j2], a1temp[j2], a2temp[j2], al1[j3], al2[j3]);exit(1);}

allcentres[k][j3]=ctemp[j2];allfactors[k][j3]=mtemp[j2]*pow(wtemp[j2],0.5);
if(a1temp[j2]!=al1[j3])
{allcentres[k][j3]=2-allcentres[k][j3];allfactors[k][j3]=-allfactors[k][j3];}

if(allfactors[k][j3]!=0){usedpreds[j3]++;}
wsums[k]+=wtemp[j2];
}

for(j=0;j<count;j++){free(ptemp[j]);}free(ptemp);
free(ctemp);free(mtemp);free(wtemp);free(a1temp);free(a2temp);
free(kindex);free(kindexb);
}	//end of k loop


if(adjust>0)	//must subtract for region predictors ???
{

}

/////////

if(num_regs>0)	//now read in regions
{
printf("Reading details for regions from %s\n", regfile);
count=countrows(regfile)-1;

ptemp=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){ptemp[j]=malloc(sizeof(char)*100);}
a1temp=malloc(sizeof(char)*count);a2temp=malloc(sizeof(char)*count);
ctemp=malloc(sizeof(float)*count);ftemp=malloc(sizeof(float)*count*num_regs);

kindex=malloc(sizeof(int)*count);kindexb=malloc(sizeof(int)*count);

//open, check then skip header line
if((input2=fopen(regfile,"r"))==NULL)
{printf("Error opening %s\n\n", regfile);exit(1);}

if(fscanf(input2, "%s ", readstring)!=1)
{printf("Error reading first element of %s\n\n", regfile);exit(1);}
if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s)\n\n", regfile, readstring);exit(1);}
readchar=0;
while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input2, "%s %c %c %f ", ptemp[j], a1temp+j, a2temp+j, ctemp+j)!=4)
{printf("Error reading first 4 elements of Row %d of %s\n\n", j+2, regfile);exit(1);}
readint--;

if(a1temp[j]==a2temp[j])
{printf("Error reading %s, both alleles for %s are the same (%c)\n\n", regfile, ptemp[j], a1temp[j]);exit(1);}

for(r=0;r<num_regs;r++)
{
if(fscanf(input2, "%f ", ftemp+j+r*count)!=1)
{printf("Error reading Region %d effect from Row %d of %s\n", r+1, j+2, regfile);exit(1);}
}
}

count2=find_strings(ptemp, count, prednames, num_preds_use, kindex, kindexb);
if(count2<count)
{printf("Data are available for %d of the %d predictors\n", count2, count);}
fclose(input2);

//now check valid and load up
for(j=0;j<count2;j++)
{
j2=kindex[j];j3=kindexb[j];
if(a1temp[j2]!=al1[j3]&&a1temp[j2]!=al2[j3])
{printf("Error reading %s, alleles for %s (%c and %c) not consistent with those in data (%c and %c)\n\n", regfile, ptemp[j2], a1temp[j2], a2temp[j2], al1[j3], al2[j3]);exit(1);}
if(a2temp[j2]!=al1[j3]&&a2temp[j2]!=al2[j3])
{printf("Error reading %s, alleles for %s (%c and %c) not consistent with those in data (%c and %c)\n\n", regfile, ptemp[j2], a1temp[j2], a2temp[j2], al1[j3], al2[j3]);exit(1);}

for(r=0;r<num_regs;r++)
{
allcentres[num_kins+r][j3]=ctemp[j2];allfactors[num_kins+r][j3]=ftemp[j2+r*count];
if(a1temp[j2]!=al1[j3])
{allcentres[num_kins+r][j3]=2-allcentres[num_kins+r][j3];allfactors[num_kins+r][j3]=-allfactors[num_kins+r][j3];}
}
usedpreds[j3]++;
}

for(j=0;j<count;j++){free(ptemp[j]);}free(ptemp);
free(ctemp);free(ftemp);free(a1temp);free(a2temp);
free(kindex);free(kindexb);
}

/////////

//find out which have been used
count=0;
for(j=0;j<num_preds_use;j++)
{
if(usedpreds[j]>0)	//used
{
keeppreds_use[count]=keeppreds[j];
for(k=0;k<num_kins+num_regs;k++)
{allcentres[k][count]=allcentres[k][j];allfactors[k][count]=allfactors[k][j];}
count++;
}
}

free(usedpreds);
printf("\n");

return(count);
}	//end of set_up_blup1

///////////////////////////

int set_up_blup2(double **mG, double **mG2, char *blupfile, char **ids1, char **ids2, int num_samples_use, int num_kins, int num_regs)
{
int i, k, r, count;
int found, *kindex;

char **ids1b, **ids2b;

float readfloat, *gtemp, *gtemp2;

FILE *input;

count=countrows(blupfile);
ids1b=malloc(sizeof(char*)*count);
ids2b=malloc(sizeof(char*)*count);
for(i=0;i<count;i++)
{ids1b[i]=malloc(sizeof(char)*500);ids2b[i]=malloc(sizeof(char)*500);}

gtemp=malloc(sizeof(float)*count*num_kins);
gtemp2=malloc(sizeof(float)*count*num_kins);

kindex=malloc(sizeof(int)*count);

if((input=fopen(blupfile,"r"))==NULL)
{printf("Error opening file %s\n\n",blupfile);exit(1);}

for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", ids1b[i], ids2b[i])!=2)
{printf("Error reading the ids on Row %d of %s\n\n", i+1, blupfile);exit(1);}

for(k=0;k<num_kins;k++)
{
if(fscanf(input, "%f %f ", gtemp2+i+k*count, gtemp+i+k*count)!=2)
{printf("Error reading random effects %d for Individual %s %s in %s\n\n", k+1, ids1b[i], ids2b[i], blupfile);exit(1);}
}
for(r=0;r<num_regs;r++)	//skip to end of row
{
if(fscanf(input, "%f %f ", &readfloat, &readfloat)!=2)
{printf("Error reading random effects %d for Individual %s %s in %s\n\n", num_kins+r+1, ids1b[i], ids2b[i], blupfile);exit(1);}
}
}
fclose(input);

//should be able to find all these individuals in famfile
found=find_ids(kindex, ids1b, ids2b, count, ids1, ids2, num_samples_use);
if(found>0){printf("Error, Ind %s %s is in %s but not in data\n\n", ids1[found-1], ids2[found-1], blupfile);exit(1);}
if(found<0){printf("Error, Ind %s %s features in data more than once\n\n", ids1[-found-1], ids2[-found-1]);exit(1);}

//fill up, starting by setting all to zero
for(i=0;i<num_samples_use;i++)
{
for(k=0;k<num_kins;k++){mG[k][i]=0;mG2[k][i]=0;}
}

for(i=0;i<count;i++)
{
for(k=0;k<num_kins;k++)
{mG[k][kindex[i]]=gtemp[i+k*count];
mG2[k][kindex[i]]=gtemp2[i+k*count];}
}

for(i=0;i<count;i++){free(ids1b[i]);free(ids2b[i]);}free(ids1b);free(ids2b);
free(gtemp);free(gtemp2);free(kindex);

return(0);
}	//end of set_up_blup2

///////////////////////////

int write_blups(double **effects, int length, double **preds, int num_samples_use, int num_kins, int num_regs, char **prednames, char *al1, char *al2, char **ids1, char **ids2, float **allcentres, char *outfile)
{
int i, j, k, r;

FILE *output1, *output2;
char filename1[500], filename2[500];


sprintf(filename1,"%s.blup",outfile);
if((output1=fopen(filename1,"w"))==NULL)
{printf("Error opening %s\n",filename1);exit(1);}
fprintf(output1, "Predictor\tA1\tA2\tTotal\t");
for(k=0;k<num_kins;k++){fprintf(output1,"K%d_Centre\tK%d_Effect\t", k+1, k+1);}
for(r=0;r<num_regs;r++){fprintf(output1,"R%d_Centre\tR%d_Effect\t", r+1, r+1);}
fprintf(output1, "\n");
for(j=0;j<length;j++)
{
if(effects[num_kins+num_regs][j]!=0)
{
fprintf(output1, "%s\t%c\t%c\t%f\t", prednames[j], al1[j], al2[j], effects[num_kins+num_regs][j]);
for(k=0;k<num_kins+num_regs;k++){fprintf(output1, "%f\t%f\t", allcentres[k][j], effects[k][j]);}
fprintf(output1, "\n");
}
}
fclose(output1);

sprintf(filename2,"%s.pred",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n",filename2);exit(1);}
fprintf(output2,"ID1\tID2\tTotal\t");
for(k=0;k<num_kins;k++){fprintf(output2,"K%d_Cont\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output2,"R%d_Cont\t", r+1);}
fprintf(output2,"\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output2, "%s\t%s\t%f\t", ids1[i], ids2[i], preds[num_kins+num_regs][i]);
for(k=0;k<num_kins+num_regs;k++){fprintf(output2, "%f\t", preds[k][i]);}
fprintf(output2, "\n");
}
fclose(output2);

return(0);
}	//end of write_blups

///////////////////////////

int write_scores(double **preds, int num_samples_use, int num_regs, char **ids1, char **ids2, char *outfile)
{
int i, r;

FILE *output;
char filename[500];


sprintf(filename,"%s.profile",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}
fprintf(output,"ID1\tID2\t");
for(r=0;r<num_regs;r++){fprintf(output,"Profile%d\t", r+1);}
fprintf(output,"\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
for(r=0;r<num_regs;r++){fprintf(output, "%f\t", preds[r][i]);}
fprintf(output, "\n");
}
fclose(output);

return(0);
}	//end of write_scores


