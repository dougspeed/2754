/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int set_up_partitions(char *folder, int num_preds_use, int partition_length, char *partpref, int num_parts, int *chr, int bychr, char **prednames, char *datafile, char *mapfile, char *bpredfile)
{
int j, count;
int misscount, start, end;
int partition, length_partitions;
int *pindex, *usedpreds;

char filename[500], partfile[500];
FILE *output;

if(strcmp(partpref,"blank")==0)	//get num of partitions and lengths (but won't use if bychr=1)
{
num_parts=(num_preds_use-1)/partition_length+1;
length_partitions=(num_preds_use-1)/num_parts+1;
}
else		//read and check the partition lists
{
pindex=malloc(sizeof(int)*(1+num_preds_use));
usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

length_partitions=0;
for(partition=0;partition<num_parts;partition++)
{
sprintf(partfile, "%s%d", partpref, partition+1);
count=read_partition(partfile, pindex, prednames, num_preds_use, partition+1);
for(j=0;j<count;j++){usedpreds[pindex[j]]++;}
if(count>length_partitions){length_partitions=count;}
}

count=0;misscount=0;
for(j=0;j<num_preds_use;j++)
{
if(usedpreds[j]>0){count++;}
if(usedpreds[j]>1)
{
if(misscount<10)
{printf("Warning, Predictor %s specified in more than one partition list\n", prednames[j]);}
misscount++;
}
}

if(misscount>0)
{printf("In total, %d predictors were specified in more than one partition list.\nThese lists are usually (but not necessarily) non-overlapping.\n\n", misscount);}

if(count==num_preds_use)
{printf("The partition lists contained all %d predictors\n\n", count);}
else
{printf("The partition lists contained only %d out of %d predictors\n\n", count, num_preds_use);}

free(pindex);free(usedpreds);
}

//open file which stores partition details and write header lines
sprintf(filename,"%spartition_details.txt",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to file \"%s\" - check write permission granted\n\n",filename);exit(1);}

fprintf(output, "Datafiles: %s %s\n", datafile, mapfile);
if(strcmp(bpredfile,"blank")==0){fprintf(output, "Using All Predictors\n");}
else{fprintf(output, "Extracting from %s\n", bpredfile);}
fprintf(output,"Partition Start_Predictor End_Predictor\n");

if(strcmp(partpref,"blank")==0)
{
if(bychr==0)	//divide into regions equally
{
for(partition=0;partition<num_parts;partition++)
{
//partitions will contain the boundaries of 1:length_partitions predictors
start=length_partitions*partition;
end=length_partitions*(partition+1);
if(end>num_preds_use){end=num_preds_use;}
fprintf(output,"%d %d %d\n", partition+1, start+1, end);
}	//end of partitions loop
}
else	//divide by chromosome
{
length_partitions=0;
j=0;count=0;num_parts=0;
while(1)
{
//should always be a new chromsome
count++;

//see if missed any chromosomes
while(count<chr[j])
{printf("No predictors on Chromosome %d\n", count);count++;}

//now find start and end
start=j;end=j;
while(chr[j]==count)
{
j++;end++;
if(j==num_preds_use){break;}
}
if(end-start>length_partitions){length_partitions=end-start;}
fprintf(output, "%d %d %d\n", num_parts+1, start+1, end);
printf("%d predictors on Chromosome %d (correspond to %skinship%d)\n", end-start, chr[start], folder, num_parts+1);
num_parts++;

if(j==num_preds_use){break;}
}	//end of while loop

}	//end of dividing by chromosome
}	//end of partitions not supplied

else	//divide predictors according to partition lists
{
for(partition=0;partition<num_parts;partition++)
{fprintf(output, "%d -1 %s%d\n", partition+1, partpref, partition+1);}
}

fclose(output);

if(strcmp(partpref,"blank")==0)
{
if(bychr==0)
{printf("The %d predictors to be examined have been split into %d partitions of (approx) length %d\n", num_preds_use, num_parts, length_partitions);}
else
{printf("\nThe %d predictors to be examined have been split into %d partitions according to chromosome, the largest of which contains %d predictors\n", num_preds_use, num_parts, length_partitions);}
}
else
{
printf("The %d predictors to be examined have been split into %d partitions according to lists %s1, ..., %s%d, the largest of which contains %d predictors\n", num_preds_use, num_parts, partpref, partpref, num_parts, length_partitions);
}

printf("The details for partitions are stored in %spartition_details.txt\n\n", folder);

return(0);
}	//end of set_up_partitions

////////////////////////////////

int get_partition_boundaries(char *folder, int partition, int *keeppreds_use, int *keeppreds, int num_preds_use, char *mapfile)
{
int j, count;
int start, end;
char **allsnps;

int readint;
char readchar;

char filename[500], partfile[500];
FILE *input;


//first set the details filename and check partition is reasonable
sprintf(filename,"%spartition_details.txt", folder);
if(partition>countrows(filename)-3)
{printf("The partition number specified (%d) is larger than the total number of partitions (%d) listed in %s\n\n", partition, countrows(filename)-2, filename);exit(1);}

//open the file and skip first partition+2 rows
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
count=0;
while(count<partition+2)
{
readchar=0;while(readchar!=10){(void)fscanf(input,"%c", &readchar);}
count++;
}

//read the row
if(fscanf(input, "%d %d ", &readint, &start)!=2)
{printf("Error reading first two details for Partition %d (Row %d) from %s\n\n", partition, partition+3, filename);exit(1);}
if(readint!=partition){printf("Error reading %s; Row %d should contain details for partition %d (not %d)\n\n", filename, partition+3, partition, readint);exit(1);}

if(start!=-1)	//start and end provided
{
start=start-1;
if(fscanf(input, "%d ", &end)!=1)
{printf("Error reading end predictor for Partition %d (row %d) from %s\n\n", partition, partition+3, filename);exit(1);}

//check sensible
if(end>num_preds_use)
{printf("Error reading %s; the end predictor (%d) of Partition %d is greater than the total number of predictors being considered (%d). Has file been altered since creation using argument \"--cut-kins\"\n\n", filename, end, partition, num_preds_use);exit(1);}

//load up keeppreds_use
for(j=0;j<end-start;j++){keeppreds_use[j]=keeppreds[start+j];}
}
else	//partition defined by the list file
{
if(fscanf(input, "%s ", partfile)!=1)
{printf("Error reading filename for Partition %d (Row %d) from %s\n\n", partition, partition+3, filename);exit(1);}

//read the list then put into keeppreds_use
count=countrows(mapfile);
allsnps=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){allsnps[j]=malloc(sizeof(char)*100);}
read_strings(mapfile, allsnps, count, 2, 0);
start=0;
end=read_partition(partfile, keeppreds_use, allsnps, count, partition);
}

fclose(input);

return(end-start);
}	//end of get_partition_boundaries

////////////////////////////////

int write_kins(char *outfile, double *kins, double *divs, float weightsum, char **ids1, char **ids2, int num_samples_use, int *keeppreds_use, char **prednames, double *centres, double *mults, float *weights, char *al1, char *al2, int data_length, int kinraw, int kingz)
{
int i, i2, j, count;
float sum, trace, value, value2;

char cmd[500];
char filename1[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500];
FILE *output1, *output2, *output3, *output4, *output5, *output6, *output7;


//set filenames
sprintf(filename1, "%s.grm.id", outfile);
sprintf(filename2, "%s.grm.bin", outfile);
sprintf(filename3, "%s.grm.N.bin", outfile);
sprintf(filename4, "%s.grm.details", outfile);
sprintf(filename5, "%s.grm.adjust", outfile);
sprintf(filename6, "%s.grm.raw", outfile);
sprintf(filename7, "%s.grm", outfile);

//write ids, kins and dividers

if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename2);exit(1);}
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename3);exit(1);}

sum=0;trace=0;
value2=weightsum;
for(i=0;i<num_samples_use;i++)

{
for(i2=0;i2<=i;i2++)
{
if(divs!=NULL){value2=divs[i+i2*num_samples_use];}
value=kins[i+i2*num_samples_use]/value2;
fwrite(&value, sizeof(float), 1, output2);
fwrite(&value2, sizeof(float), 1, output3);

if(i2<i){sum+=2*value;}
else{sum+=value;trace+=value;}
}
}
fclose(output2);
fclose(output3);

if((output1=fopen(filename1,"w"))==NULL)
{printf("Error writing to %s\n\n",filename1);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output1, "%s %s\n", ids1[i], ids2[i]);}
fclose(output1);

if(keeppreds_use!=NULL)	//write details
{
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s\n\n",filename4);exit(1);}
fprintf(output4, "Predictor Index Centre Scaling Weight A1 A2\n");
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-1&&weights[j]>0)
{fprintf(output4, "%s %d %f %f %f %c %c\n", prednames[j], keeppreds_use[j]+1, centres[j], mults[j], weights[j], al1[j], al2[j]);
count++;}
}
fclose(output4);
}

//and adjustments
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s\n\n",filename5);exit(1);}
if(keeppreds_use!=NULL)
{fprintf(output5, "Mean_Trace: %.3f\nMean_Kinship: %.3f\nNum_Preds %d\nSum_Weights %.2f\n", trace/num_samples_use, sum/num_samples_use/num_samples_use, count, weightsum);}
else
{fprintf(output5, "Mean_Trace: %.3f\nMean_Kinship: %.3f\nNum_Preds NA\nSum_Weights %.2f\n", trace/num_samples_use, sum/num_samples_use/num_samples_use, weightsum);}
fclose(output5);

if(kinraw==1)
{
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s\n\n",filename6);exit(1);}

for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<num_samples_use;i2++)
{
value2=weightsum;
if(divs!=NULL){value2=divs[i+i2*num_samples_use];}
fprintf(output6, "%f ", kins[i+i2*num_samples_use]/value2);}
fprintf(output6, "\n");
}
fclose(output6);
}

if(kingz==1)
{
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s\n\n",filename7);exit(1);}

value=weightsum;
value2=value+0.5;
if(value2<1){value2=1;}
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<=i;i2++)
{
if(divs!=NULL)
{
value=divs[i+i2*num_samples_use];
value2=value+0.5;
if(value2<1){value2=1;}
}
fprintf(output7, "%d %d %d %f\n", i+1, i2+1, (int)value2, kins[i+i2*num_samples_use]/value);}
}
fclose(output7);

//try to compress the kinship file
sprintf(cmd, "gzip -f %s", filename7);
system(cmd);
}

return(0);
}	//end of write_kins

////////////////////////////////

float adjust_details(char *kinstem, char **prednames, double *centres, double *mults, float *weights, char *al1, char *al2, int *keeppreds_use, int rdata_length, char **rprednames, double *rcentres, double *rmults, float *rweights, char *ral1, char *ral2)
{
int j, j2, j3, count, count2;
int *kindex, *kindexb;
float weightsum;

int readint;
float readcentre, readmult;
char readchar;

char filename[500];
FILE *input;

sprintf(filename, "%s.grm.details", kinstem);
count=countrows(filename)-1;

//open and skip header line
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input,"%c", &readchar);}

weightsum=0;
for(j=0;j<count;j++)
{
if(fscanf(input, "%s %d %f %f %f %c %c ", prednames[j], &readint, &readcentre, &readmult, weights+j, al1+j, al2+j)!=7)
{printf("Error reading Row %d of %s\n\n", j+2, filename);exit(1);}
keeppreds_use[j]=readint-1;
centres[j]=readcentre;
mults[j]=readmult;
weightsum+=weights[j];

if(al1[j]==al2[j])
{printf("Error reading %s; both alleles for %s are the same (%c)\n\n", filename, prednames[j], al1[j]);exit(1);}
}
fclose(input);

//now correct the region details
kindex=malloc(sizeof(int)*rdata_length);
kindexb=malloc(sizeof(int)*rdata_length);
count2=find_strings(rprednames, rdata_length, prednames, count, kindex, kindexb);
if(count2==0){printf("Error, none of the %d region predictors feature in %s\n\n", rdata_length, filename);exit(1);}
if(count2==count){printf("Error, there will be no predictors left after removing the %d region predictors found in %s\n\n", count, filename);exit(1);}
printf("%d of the %d region predictors are present in %s; these will be removed\n\n", count2, rdata_length, filename);

for(j=0;j<count2;j++)
{
j2=kindex[j];j3=kindexb[j];
rcentres[j2]=centres[j3];
rmults[j2]=mults[j3];
rweights[j2]=weights[j3];
weights[j3]=0;
}

free(kindex);free(kindexb);

return(weightsum);
}	//end of adjust_details

////////////////////////////////

int manip_kins(char *outfile, int num_kins, char **kinstems, int samedata, int kindetails, int kinraw, int kingz, int type)
//type=0, join-kins; type=1, add; type=2, subtract
{
int i, j, j2, j3, k, count, count2, count3;
int num_samples_use, num_preds, num_preds_approx, *kindex, *kindexb, *pindex, *usedpreds, *founds;
float weightsum;

int readint;
char readchar, readstring[100];

char **ids1, **ids2, **ids1b, **ids2b, **usesamps, **allsamps, **allsnps, **allsnpsb;

int *indexes, *itemp;
float *ctemp, *mtemp, *wtemp, *weights;
double *centres, *mults;
char **ptemp, *a1temp, *a2temp, *al1, *al2;

double *kins, *divs;

char filename1[500], filename2[500];
FILE *input2;


//first work out which ids we will be using

if(type==0)	//all ids should match, so get from the first file and check remainder
{
sprintf(filename1, "%s.grm.id", kinstems[0]);
num_samples_use=countrows(filename1);
ids1=malloc(sizeof(char*)*num_samples_use);ids2=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{ids1[i]=malloc(sizeof(char)*100);ids2[i]=malloc(sizeof(char)*100);}
read_famfile(filename1, ids1, ids2, num_samples_use, NULL);

for(k=1;k<num_kins;k++)
{
sprintf(filename1, "%s.grm.id", kinstems[k]);
count=countrows(filename1);
if(count!=num_samples_use)
{printf("Error, number of ids in %s (%d) does not match number in %s.grm.id (%d)\nIt seems different keepfiles were used for different partitions\n\n", filename1, count, kinstems[0], num_samples_use);exit(1);}

ids1b=malloc(sizeof(char*)*count);ids2b=malloc(sizeof(char*)*count);
for(i=0;i<count;i++)
{ids1b[i]=malloc(sizeof(char)*100);ids2b[i]=malloc(sizeof(char)*100);}
read_famfile(filename1, ids1b, ids2b, count, NULL);

for(i=0;i<num_samples_use;i++)
{
if(strcmp(ids1[i],ids1b[i])!=0||strcmp(ids2[i],ids2b[i])!=0)
{printf("Error, ids in %s do not match those in %s.grm.id\n\n", filename1, kinstems[0]);}
}
for(i=0;i<count;i++){free(ids1b[i]);free(ids2b[i]);}free(ids1b);free(ids2b);
}
}

if(type==1||type==2)	//find overlap of ids
{
//start with first kinship file
sprintf(filename1, "%s.grm.id", kinstems[0]);
count=countrows(filename1);
allsamps=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){allsamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(filename1, allsamps, count);

founds=malloc(sizeof(int)*count);
for(i=0;i<count;i++){founds[i]=1;}
kindex=malloc(sizeof(int)*count);

for(k=1;k<num_kins;k++)	//which are found in each kinship file
{
sprintf(filename1, "%s.grm.id", kinstems[k]);
count2=countrows(filename1);
usesamps=malloc(sizeof(char*)*count2);
for(i=0;i<count2;i++){usesamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(filename1, usesamps, count2);

count3=find_strings(allsamps, count, usesamps, count2, kindex, NULL);
for(i=0;i<count3;i++){founds[kindex[i]]++;}
for(i=0;i<count2;i++){free(usesamps[i]);}free(usesamps);
}

for(i=0;i<count;i++){free(allsamps[i]);}free(allsamps);

num_samples_use=0;
for(i=0;i<count;i++)
{
if(founds[i]==num_kins){kindex[num_samples_use]=i;num_samples_use++;}
}
printf("%d ids are common to the %d kinship matrices\n\n", num_samples_use, num_kins);
ids1=malloc(sizeof(char*)*num_samples_use);ids2=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{ids1[i]=malloc(sizeof(char)*100);ids2[i]=malloc(sizeof(char)*100);}

sprintf(filename1, "%s.grm.id", kinstems[0]);
read_famfile(filename1, ids1, ids2, num_samples_use, kindex);
free(founds);free(kindex);
}

///////////////////////////

if(kindetails==1)	//now sort details and get weightsum
{
printf("Sorting details\n");

if(samedata==1)	//then indexes should be consistent across detail files
{
//first set num_preds to max index observed
num_preds=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename2, "%s.grm.details", kinstems[k]);
count=countrows(filename2)-1;
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input2, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row %d of %s\n\n", j+2, filename2);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}
if(readint<=0)
{printf("Error reading %s; Index in Row %d is not positive (%d)\n\n", filename2, j+2, readint);exit(1);}
if(readint>num_preds){num_preds=readint;}
}
fclose(input2);
}

//now fill names of predictors used
allsnps=malloc(sizeof(char*)*num_preds);
for(j=0;j<num_preds;j++){allsnps[j]=malloc(sizeof(char)*100);strcpy(allsnps[j],"blank");}
for(k=0;k<num_kins;k++)
{
sprintf(filename2, "%s.grm.details", kinstems[k]);
count=countrows(filename2)-1;
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input2, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row %d of %s\n\n", j+2, filename2);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}
readint--;

if(strcmp(allsnps[readint],readstring)!=0&&strcmp(allsnps[readint],"blank")!=0)
{printf("Error reading %s; Index %d refers to both Predictor %s and %s\n\n", filename2, readint+1, allsnps[readint], readstring);exit(1);}
strcpy(allsnps[readint],readstring);
}
fclose(input2);
}
}	//end of samedata=1

if(samedata==0)	//then get list of unique predictors
{
//start by reading in list of all predictors (ignoring whether duplicates)
num_preds_approx=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename2, "%s.grm.details", kinstems[k]);
count=countrows(filename2)-1;
num_preds_approx+=count;
}
allsnpsb=malloc(sizeof(char*)*num_preds_approx);
for(j=0;j<num_preds_approx;j++){allsnpsb[j]=malloc(sizeof(char)*100);}
num_preds_approx=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename2, "%s.grm.details", kinstems[k]);
count=countrows(filename2)-1;
read_strings(filename2, allsnpsb+num_preds_approx, count, 1, 1);
num_preds_approx+=count;
}

//sort this list and find number of unique
pindex=malloc(sizeof(int)*num_preds_approx);
sortindex_string(allsnpsb, pindex, num_preds_approx);
num_preds=1;
for(j=1;j<num_preds_approx;j++)
{
if(strcmp(allsnpsb[pindex[j-1]],allsnpsb[pindex[j]])!=0)
num_preds++;
}
allsnps=malloc(sizeof(char*)*num_preds);
for(j=0;j<num_preds;j++){allsnps[j]=malloc(sizeof(char)*100);}
strcpy(allsnps[0],allsnpsb[pindex[0]]);num_preds=1;
for(j=1;j<num_preds_approx;j++)
{
if(strcmp(allsnpsb[pindex[j-1]],allsnpsb[pindex[j]])!=0)
{strcpy(allsnps[num_preds],allsnpsb[pindex[j]]);num_preds++;}
}

for(j=0;j<num_preds_approx;j++){free(allsnpsb[j]);}free(allsnpsb);
free(pindex);
}	//end of samedata=0

/////////

//allocate then read in
indexes=malloc(sizeof(int)*num_preds);
centres=malloc(sizeof(double)*num_preds);
mults=malloc(sizeof(double)*num_preds);
weights=malloc(sizeof(float)*num_preds);
al1=malloc(sizeof(char)*num_preds);
al2=malloc(sizeof(char)*num_preds);

for(j=0;j<num_preds;j++){weights[j]=0;}

usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-1;}

for(k=0;k<num_kins;k++)
{
sprintf(filename2, "%s.grm.details", kinstems[k]);
count=countrows(filename2)-1;

ptemp=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){ptemp[j]=malloc(sizeof(char)*100);}
ctemp=malloc(sizeof(float)*count);itemp=malloc(sizeof(int)*count);
mtemp=malloc(sizeof(float)*count);wtemp=malloc(sizeof(float)*count);
a1temp=malloc(sizeof(char)*count);a2temp=malloc(sizeof(char)*count);

kindex=malloc(sizeof(int)*count);kindexb=malloc(sizeof(int)*count);

//open and skip header line
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input2,"%c", &readchar);}

count2=0;
for(j=0;j<count;j++)
{
if(fscanf(input2, "%s %d %f %f %f %c %c ", ptemp[j], itemp+j, ctemp+j, mtemp+j, wtemp+j, a1temp+j, a2temp+j)!=7)
{printf("Error reading Row %d of %s\n\n", j+2, filename2);exit(1);}
readint=itemp[j]-1;

if(a1temp[j]==a2temp[j])
{printf("Error reading %s; both alleles for %s are the same (%c)\n\n", filename2, ptemp[j], a1temp[j]);exit(1);}

if(samedata==1)	//then should be able to use indexes
{
if(strcmp(ptemp[j], allsnps[readint])!=0){printf("Doug error 3\n");exit(1);}
kindex[count2]=j;kindexb[count2]=readint;count2++;
}
}
fclose(input2);

if(samedata==0)	//can not use indexes, so find manually
{
count2=find_strings(ptemp, count, allsnps, num_preds, kindex, kindexb);
if(count2<count)
{printf("Doug Error %d %d\n", count2, count);}
}

//now load up
for(j=0;j<count;j++)
{
j2=kindex[j];j3=kindexb[j];

if(usedpreds[j3]==-1)	//it's new
{
if(type==2&&k!=0)
{printf("Error reading %s; can not subtract Predictor %s as it was not in %s.grm.details\n\n", filename2, allsnps[j3], kinstems[0]);exit(1);}

centres[j3]=ctemp[j2];indexes[j3]=itemp[j2]-1;
mults[j3]=mtemp[j2];weights[j3]=wtemp[j2];
al1[j3]=a1temp[j2];al2[j3]=a2temp[j2];
usedpreds[j3]=k;
}
else	//already used (so can't be looking at first)
{
if(a1temp[j2]!=al1[j3]&&a1temp[j2]!=al2[j3])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename2, ptemp[j2], a1temp[j2], a2temp[j2], kinstems[usedpreds[j3]], al1[j3], al2[j3]);exit(1);}
if(a2temp[j2]!=al1[j3]&&a2temp[j2]!=al2[j3])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename2, ptemp[j2], a1temp[j2], a2temp[j2], kinstems[usedpreds[j3]], al1[j3], al2[j3]);exit(1);}

if(a1temp[j2]!=al1[j3])	//flip
{ctemp[j2]=2-ctemp[j2];}
if(ctemp[j2]!=centres[j3]||mtemp[j2]!=mults[j3])	//should match
{printf("Error reading %s; centres and scaling for %s different to those in %s.grm.details\nIf a predictor features in more than one kinship file, the same individuals must be used when constructing each\n\n", filename2, allsnps[j3], kinstems[usedpreds[j3]]);exit(1);}

if(itemp[j2]!=indexes[j3]+1&&samedata==1)
{printf("Doug error 5\n\n");exit(1);}

if(type==0||type==1)
{weights[j3]+=wtemp[j2];}
else
{
weights[j3]-=wtemp[j2];
if(weights[j3]<0)
{printf("Error reading %s, predictor %s has been subtracted \"too much\"\n\n", filename2, allsnps[j3]);exit(1);}
}
}	//end of already used
}	//end of j loop

for(j=0;j<count;j++){free(ptemp[j]);}free(ptemp);
free(ctemp);free(itemp);free(mtemp);free(wtemp);free(a1temp);free(a2temp);
free(kindex);free(kindexb);
}	//end of k loop

if(samedata==0)	//indexes can not be used so set all to -2
{
for(j=0;j<num_preds;j++){indexes[j]=-2;}
}
}	//end of kindetails=1

////////////////////////////////

printf("Sorting kinships\n");
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
divs=malloc(sizeof(double)*num_samples_use*num_samples_use);
for(i=0;i<num_samples_use*num_samples_use;i++){kins[i]=0;divs[i]=0;}

for(k=0;k<num_kins;k++)
{
if(k==0||type==0||type==1){read_kinfile(kinstems[k], kins, divs, num_samples_use, ids1, ids2, 1);}
else{read_kinfile(kinstems[k], kins, divs, num_samples_use, ids1, ids2, 2);}
}

//////////

if(kindetails==1)
{
weightsum=0;for(j=0;j<num_preds;j++){weightsum+=weights[j];}
write_kins(outfile, kins, divs, weightsum, ids1, ids2, num_samples_use, indexes, allsnps, centres, mults, weights, al1, al2, num_preds, kinraw, kingz);
}
else	//get a fudge as the max div recorded
{
weightsum=0;
for(i=0;i<num_samples_use*num_samples_use;i++)
{
if(divs[i]>weightsum){weightsum=divs[i];}
}
write_kins(outfile, kins, divs, weightsum, ids1, ids2, num_samples_use, NULL, NULL, NULL, NULL, NULL, NULL, NULL, num_preds, kinraw, kingz);
}

for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);}free(ids1);free(ids2);
for(j=0;j<num_preds;j++){free(allsnps[j]);}free(allsnps);
free(indexes);free(centres);free(mults);free(weights);free(al1);free(al2);
free(usedpreds);free(kins);free(divs);

return(0);
}	//end of manip_kins

////////////////////////////////









