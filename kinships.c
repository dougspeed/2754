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

int write_kins(char *outfile, double *kins, float weightsum, char **ids1, char **ids2, int num_samples_use, int *keeppreds_use, char **prednames, double *centres, double *mults, float *weights, char *al1, char *al2, int data_length, int kinraw, int kingz)
{
int i, i2, j, count;
float sum, trace, value;

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
if((output1=fopen(filename1,"w"))==NULL)
{printf("Error writing to %s\n\n",filename1);exit(1);}
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename2);exit(1);}
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename3);exit(1);}

sum=0;trace=0;
for(i=0;i<num_samples_use;i++)
{
fprintf(output1, "%s %s\n", ids1[i], ids2[i]);
for(i2=0;i2<=i;i2++)
{
value=kins[i+i2*num_samples_use]/weightsum;
fwrite(&value, sizeof(float), 1, output2);
fwrite(&weightsum, sizeof(float), 1, output3);

if(i2<i){sum+=2*value;}
else{sum+=value;trace+=value;}
}
}
fclose(output1);
fclose(output2);
fclose(output3);

//write details
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

//and adjustments
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s\n\n",filename5);exit(1);}
fprintf(output5, "Mean_Trace: %.3f\nMean_Kinship: %.3f\nNum_Preds %d\nSum_Weights %.2f\n", trace/num_samples_use, sum/num_samples_use/num_samples_use, count, weightsum);

if(kinraw==1)
{
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s\n\n",filename6);exit(1);}

for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<num_samples_use;i2++)
{fprintf(output6, "%f ", kins[i+i2*num_samples_use]/weightsum);}
fprintf(output6, "\n");
}
fclose(output6);
}

if(kingz==1)
{
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s\n\n",filename7);exit(1);}

value=weightsum+0.5;
if(value<1){value=1;}
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<=i;i2++)
{fprintf(output7, "%d %d %d %f\n", i+1, i2+1, (int)value, kins[i+i2*num_samples_use]/weightsum);}
}
fclose(output7);

//try to compress the kinship file
sprintf(cmd, "gzip -f %s", filename7);
system(cmd);
}

return(0);
}	//end of write_kins

////////////////////////////////

int manip_kins(char *outfile, int num_kins, char **kinstems, char *mapfile, int kinraw, int kingz, int type)
//type=0, add; type=1, subtract
{
int i, i2, j, j2, j3, k, count, count2, flag, flag2;
int num_samples_use, num_preds, num_preds_approx, *kindex, *kindexb, *pindex, *usedpreds;
float sum, trace, value, value2, weightsum;

int readint;
char readchar;

char **ids1, **ids2, **usesamps, **allsamps, **allsnps, **allsnpsb;

int *indexes, *itemp;
float *centres, *mults, *weights, *ctemp, *mtemp, *wtemp;
char **ptemp, *al1, *al2, *a1temp, *a2temp;

double *kins, *divs;

char cmd[500];
char filename1[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500];
FILE *output2, *input3, *output4, *output5, *output6, *output7, *output8, *output9;


//get ids from first file, check others match
sprintf(filename1, "%s.grm.id", kinstems[0]);
num_samples_use=countrows(filename1);
ids1=malloc(sizeof(char*)*num_samples_use);ids2=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{ids1[i]=malloc(sizeof(char)*100);ids2[i]=malloc(sizeof(char)*100);}
read_famfile(filename1, ids1, ids2, num_samples_use, NULL);

usesamps=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++){usesamps[i]=malloc(sizeof(char)*200);}
for(i=0;i<num_samples_use;i++){sprintf(usesamps[i], "%s___%s", ids1[i], ids2[i]);}

for(k=1;k<num_kins;k++)
{
sprintf(filename1, "%s.grm.id", kinstems[k]);
count=countrows(filename1);
if(count!=num_samples_use)
{printf("Error, number of ids in %s (%d) does not match number in %s.grm.id (%d)\n\n", filename1, count, kinstems[0], num_samples_use);exit(1);}

allsamps=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){allsamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(filename1, allsamps, count);

kindex=malloc(sizeof(int)*count);
count2=find_strings(allsamps, count, usesamps, num_samples_use, kindex, NULL);
if(count2!=num_samples_use)
{printf("Error reading %s; only found %d of the %d individuals in %s.grm.id\n\n", filename1, count2, num_samples_use, kinstems[0]);exit(1);}
}

//write ids
sprintf(filename2, "%s.grm.id", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output2, "%s %s\n", ids1[i], ids2[i]);}
fclose(output2);

//////////

//now sort details
if(strcmp(mapfile,"blank")!=0)	//read predictor names from mapfile
{
num_preds=countrows(mapfile);
allsnps=malloc(sizeof(char*)*num_preds);
for(j=0;j<num_preds;j++){allsnps[j]=malloc(sizeof(char)*100);}
read_strings(mapfile, allsnps, num_preds, 2, 0);
}
else	//then get list of unique predictors
{
num_preds_approx=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename3, "%s.grm.details", kinstems[k]);
count=countrows(filename3)-1;
num_preds_approx+=count;
}

allsnpsb=malloc(sizeof(char*)*num_preds_approx);
for(j=0;j<num_preds_approx;j++){allsnpsb[j]=malloc(sizeof(char)*100);}

num_preds_approx=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename3, "%s.grm.details", kinstems[k]);
count=countrows(filename3)-1;
read_strings(filename3, allsnpsb+num_preds_approx, count, 1, 1);
num_preds_approx+=count;
}

//sort and get unique list 
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
}

//allocate then read in
indexes=malloc(sizeof(int)*num_preds);
centres=malloc(sizeof(float)*num_preds);
mults=malloc(sizeof(float)*num_preds);
weights=malloc(sizeof(float)*num_preds);
al1=malloc(sizeof(char)*num_preds);
al2=malloc(sizeof(char)*num_preds);

usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-1;}

flag2=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename3, "%s.grm.details", kinstems[k]);
count=countrows(filename3)-1;

ptemp=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){ptemp[j]=malloc(sizeof(char)*100);}
ctemp=malloc(sizeof(float)*count);itemp=malloc(sizeof(int)*count);
mtemp=malloc(sizeof(float)*count);wtemp=malloc(sizeof(float)*count);
a1temp=malloc(sizeof(char)*count);a2temp=malloc(sizeof(char)*count);

kindex=malloc(sizeof(int)*count);kindexb=malloc(sizeof(int)*count);

//open and skip header line
if((input3=fopen(filename3,"r"))==NULL)
{printf("Error opening %s\n\n", filename3);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input3,"%c", &readchar);}

flag=0;count2=0;
for(j=0;j<count;j++)
{
if(fscanf(input3, "%s %d %f %f %f %c %c ", ptemp[j], itemp+j, ctemp+j, mtemp+j, wtemp+j, a1temp+j, a2temp+j)!=7)
{printf("Error reading Row %d of %s\n\n", j+2, filename3);exit(1);}
readint=itemp[j]-1;

if(a1temp[j]==a2temp[j])
{printf("Error reading %s, both alleles for %s are the same (%c)\n\n", filename3, ptemp[j], a1temp[j]);exit(1);}

if(flag==0&&readint>=0&&readint<num_preds)
{
if(strcmp(ptemp[j], allsnps[readint])==0)
{kindex[count2]=j;kindexb[count2]=readint;count2++;}
else
{flag=1;}
}
else
{flag=1;}
}
fclose(input3);

if(flag==1)	//could not find predictors using index, so find manually
{
count2=find_strings(ptemp, count, allsnps, num_preds, kindex, kindexb);
if(count2<count)
{printf("Doug Error %d %d\n", count2, count);}
}	//end of flag=1

//now load up
for(j=0;j<count;j++)
{
j2=kindex[j];j3=kindexb[j];

if(usedpreds[j3]==-1)	//it's new
{
if(type==1&&k!=0)
{printf("Error reading %s; can not subtract Predictor %s as was not in %s.grm.details\n\n", filename3, allsnps[j3], kinstems[0]);exit(1);}

centres[j3]=ctemp[j2];indexes[j3]=itemp[j2];
mults[j3]=mtemp[j2];weights[j3]=wtemp[j2];
al1[j3]=a1temp[j2];al2[j3]=a2temp[j2];
usedpreds[j3]=k;
}
else	//already used (so can't be looking at first)
{
if(a1temp[j2]!=al1[j3]&&a1temp[j2]!=al2[j3])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename3, ptemp[j2], a1temp[j2], a2temp[j2], kinstems[usedpreds[j3]], al1[j3], al2[j3]);exit(1);}
if(a2temp[j2]!=al1[j3]&&a2temp[j2]!=al2[j3])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename3, ptemp[j2], a1temp[j2], a2temp[j2], kinstems[usedpreds[j3]], al1[j3], al2[j3]);exit(1);}

if(a1temp[j2]!=al1[j3])	//flip
{ctemp[j2]=2-ctemp[j2];}
if(ctemp[j2]!=centres[j3]||mtemp[j2]!=mults[j3])	//should match
{printf("Error reading %s; centres and scaling for %s different to those in %s.grm.details\n\n", filename3, allsnps[j3], kinstems[usedpreds[j3]]);exit(1);}
if(itemp[j2]!=indexes[j3]){flag2=1;}

if(type==0)
{weights[j3]+=wtemp[j2];}
else
{
weights[j3]-=wtemp[j2];
if(weights[j3]<0)
{printf("Error reading %s, predictor %s has been subtracted \"too much\"\n\n", filename3, allsnps[j3]);exit(1);}
}
}	//end of already used
}	//end of j loop

for(j=0;j<count;j++){free(ptemp[j]);}free(ptemp);
free(ctemp);free(itemp);free(mtemp);free(wtemp);free(a1temp);free(a2temp);
free(kindex);free(kindexb);
}	//end of k loop

pindex=malloc(sizeof(int)*num_preds);
if(flag2==0)	//sort into order
{
sortindex_int(indexes, pindex, num_preds);
}
if(flag2==1)
{
printf("Warning, the kinships files were computed from more than one dataset\n\n");
for(j=0;j<num_preds;j++){indexes[j]=-1;pindex[j]=j;}
}

//write details
sprintf(filename4, "%s.grm.details", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s\n\n",filename4);exit(1);}
fprintf(output4, "Predictor Index Centre Scaling Weight A1 A2\n");
count=0;weightsum=0;
for(j=0;j<num_preds;j++)
{
j2=pindex[j];
if(weights[j2]>0)
{
fprintf(output4, "%s %d %f %f %f %c %c\n", allsnps[j2], indexes[j2], centres[j2], mults[j2], weights[j2], al1[j2], al2[j2]);
count++;
weightsum+=weights[j2];
}
}
fclose(output4);

for(j=0;j<num_preds;j++){free(allsnps[j]);}free(allsnps);
free(centres);free(indexes);free(mults);free(weights);free(al1);free(al2);
free(usedpreds);

/////////

//sort kins
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
divs=malloc(sizeof(double)*num_samples_use*num_samples_use);
for(i=0;i<num_samples_use*num_samples_use;i++){kins[i]=0;divs[i]=0;}

for(k=0;k<num_kins;k++)
{
if(k==0||type==0){read_kinfile(kinstems[k], kins, divs, num_samples_use, ids1, ids2, 1);}
else{read_kinfile(kinstems[k], kins, divs, num_samples_use, ids1, ids2, 2);}
}

//and write
sprintf(filename5, "%s.grm.bin", outfile);
if((output5=fopen(filename5,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename5);exit(1);}
sprintf(filename6, "%s.grm.N.bin", outfile);
if((output6=fopen(filename6,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename6);exit(1);}

sum=0;trace=0;
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<=i;i2++)
{
value2=divs[i+i2*num_samples_use];
value=kins[i+i2*num_samples_use]/value2;
fwrite(&value, sizeof(float), 1, output5);
fwrite(&value2, sizeof(float), 1, output6);

if(i2<i){sum+=2*value;}
else{sum+=value;trace+=value;}
}
}
fclose(output5);
fclose(output6);

//and adjustments
sprintf(filename7, "%s.grm.adjust", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s\n\n",filename7);exit(1);}
fprintf(output7, "Mean_Trace: %.3f\nMean_Kinship: %.3f\nNum_Preds %d\nSum_Weights %.2f\n", trace/num_samples_use, sum/num_samples_use/num_samples_use, count, weightsum);

if(kingz==1)
{
sprintf(filename8, "%s.grm", outfile);
if((output8=fopen(filename8,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename8);exit(1);}
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<=i;i2++)
{
value2=divs[i+i2*num_samples_use]+0.5;
if(value2<1){value2=1;}
fprintf(output8, "%d %d %d %f\n", i+1, i2+1, (int)value2, kins[i+i2*num_samples_use]/divs[i+i2*num_samples_use]);
}
}
fclose(output8);

//try to compress the kinship file
sprintf(cmd, "gzip -f %s", filename8);
system(cmd);
}

if(kinraw==1)
{
sprintf(filename9, "%s.grm.raw", outfile);
if((output9=fopen(filename9,"w"))==NULL)
{printf("Error writing to %s\n\n", filename9);exit(1);}

for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<num_samples_use;i2++)
{fprintf(output9, "%f ", kins[i+i2*num_samples_use]/divs[i+i2*num_samples_use]);}
fprintf(output9, "\n");
}
fclose(output9);
}

return(0);
}	//end of manip_kins

////////////////////////////////

int adjust_details(char *outfile, char *kinstem, char **rprednames, double *rcentres, double *rmults, float *rweights, int rlength)
{
int j, j2, j3, count;
int num_preds, *pindex, *kindex, *kindexb, rlength2;

char readchar;

int *indexes;
float *centres, *mults, *weights;
char **prednames, *al1, *al2, **rprednamesb;


char filename1[500], filename2[500];
FILE *input1, *output2;


//first read in the details
sprintf(filename1, "%s.grm.details", kinstem);
num_preds=countrows(filename1)-1;

prednames=malloc(sizeof(char*)*num_preds);
for(j=0;j<num_preds;j++){prednames[j]=malloc(sizeof(char)*100);}

indexes=malloc(sizeof(int)*num_preds);
centres=malloc(sizeof(float)*num_preds);
mults=malloc(sizeof(float)*num_preds);
weights=malloc(sizeof(float)*num_preds);
al1=malloc(sizeof(char)*num_preds);
al2=malloc(sizeof(char)*num_preds);

if((input1=fopen(filename1,"r"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
readchar=0;while(readchar!=10){(void)fscanf(input1,"%c", &readchar);}
for(j=0;j<num_preds;j++)
{
if(fscanf(input1, "%s %d %f %f %f %c %c ", prednames[j], indexes+j, centres+j, mults+j, weights+j, al1+j, al2+j)!=7)
{printf("Error reading Row %d of %s\n\n", j+2, filename1);exit(1);}
}
fclose(input1);

//reduce region predictors to those which contribute
rprednamesb=malloc(sizeof(char*)*rlength);
for(j=0;j<rlength;j++){rprednamesb[j]=malloc(sizeof(char)*100);}
pindex=malloc(sizeof(int)*rlength);
rlength2=0;
for(j=0;j<rlength;j++)
{
if(rmults[j]!=-1&&rweights[j]>0)
{
strcpy(rprednamesb[rlength2],rprednames[j]);
pindex[rlength2]=j;
rlength2++;
}
}

//now find these and eliminate
kindex=malloc(sizeof(int)*rlength);kindexb=malloc(sizeof(int)*rlength);
count=find_strings(rprednames, rlength2, prednames, num_preds, kindex, kindexb);

printf("count %d re %d\n", count, rlength2);
if(count<rlength2)
{printf("Error reading %s; unable to find all the region predictors\n\n", filename1);exit(1);} 

for(j=0;j<rlength2;j++)
{
j2=pindex[kindex[j]];j3=kindexb[j];
printf("%s and %s\n", rprednames[j2], prednames[j3]);
}






}	//end of check details











