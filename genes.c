/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int set_up_genes1(char *genefile, int chunksbp, double chunks, int *chr, double *bp, float *weights, int num_preds_use)
{
int j;
int prev, numchr;
int num_genes_approx;
double *cumsum;

//work out (upper limit for) number of genes
cumsum=malloc(sizeof(double)*num_preds_use);

if(chunks==-1&&chunksbp==-1)	//reading from genefile
{
num_genes_approx=countrows(genefile);
}
if(chunks!=-1)	//using weights
{
cumsum[0]=weights[0];
numchr=1;prev=chr[0];
for(j=1;j<num_preds_use;j++)
{
cumsum[j]=cumsum[j-1]+weights[j];
if(chr[j]!=prev){numchr++;prev=chr[j];}
}
num_genes_approx=(int)(cumsum[num_preds_use-1]/chunks*2)+2*numchr;
}
if(chunksbp!=-1)	//using basepairs
{
num_genes_approx=0;
numchr=1;prev=chr[0];
for(j=1;j<num_preds_use;j++)
{
if(chr[j]!=prev||j==num_preds_use-1)	//new chr, or end element
{num_genes_approx+=(int)(bp[j-1]/chunksbp*2)+2;prev=chr[j];}
}
}

free(cumsum);

return(num_genes_approx);
}	//end of set_up_genes1

///////////////////////////

int set_up_genes2(char *folder, int *gchr, double *gbp1, double *gbp2, char **genenames, int *gstarts, int *gends, char *genefile, int chunksbp, double chunks, int *chr, double *bp, float *weights, int num_preds_use, int gene_buffer, double minweight, int overlap, int partition_length, char *pvafile, int *keeppreds_use, char *mapfile, char *datafile, char *kpredfile, int bin)
{
int j, j2, prev, prev2, last, mark, mark2, partition;
int count, count2, count3, misscount;
int num_genes_approx, *usedpreds;
float *pvalues, weightsum;
double *cumsum, minpvalue, minpvalue2;

char filename1[500], filename2[500], filename3[500];
FILE *output1, *output2, *output3;


//get start and end basepairs (weights) of each gene
cumsum=malloc(sizeof(double)*num_preds_use);

if(chunks==-1&&chunksbp==-1)	//read gene details from genefile
{
num_genes_approx=countrows(genefile);
read_genefile(genefile, genenames, gchr, gbp1, gbp2, num_genes_approx);
}

if(chunks!=-1)	//using weights
{
cumsum[0]=weights[0];
gchr[0]=chr[0];gbp1[0]=0;gbp2[0]=chunks;
num_genes_approx=1;prev=chr[0];
for(j=1;j<num_preds_use;j++)
{
cumsum[j]=cumsum[j-1]+weights[j];
if(chr[j]!=prev||cumsum[j]>gbp1[num_genes_approx-1]+chunks/2)	//then add a new gene
{
gchr[num_genes_approx]=chr[j];
gbp1[num_genes_approx]=gbp1[num_genes_approx-1]+chunks/2;
if(chr[j]!=prev){gbp1[num_genes_approx]=cumsum[j-1];}
gbp2[num_genes_approx]=gbp1[num_genes_approx]+chunks;
if(chr[j]!=prev){prev=chr[j];}
num_genes_approx++;
}
}
}

if(chunksbp!=-1)	//using basepairs - may end up with empty genes
{
num_genes_approx=0;prev=chr[0];
for(j=0;j<num_preds_use;j++)
{
if(chr[j]!=prev||j==num_preds_use-1)	//new chr, or end element
{
for(count=0;count<(int)(bp[j-1]/chunksbp*2)+2;count++)
{
gchr[num_genes_approx]=prev;
gbp1[num_genes_approx]=count*chunksbp/2+1;
gbp2[num_genes_approx]=gbp1[num_genes_approx]+chunksbp;
num_genes_approx++;
}
prev=chr[j];
}
}
}


//sort out pvalue (might not be using)
pvalues=malloc(sizeof(float)*num_preds_use);
for(j=0;j<num_preds_use;j++){pvalues[j]=2;}

if(strcmp(pvafile,"blank")!=0)
{
//read_pvafile = give warnings for nas
}


//now see which snps inside each gene

usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

mark=0;misscount=0;
for(j=0;j<num_genes_approx;j++)
{
//first get to correct chr
if(chr[mark]<gchr[j])
{
while(chr[mark]<gchr[j]&&mark<num_preds_use-1){mark++;}
}

if(chr[mark]==gchr[j])	//have found the correct chr
{
//see if we should backtrack
while(mark>0)
{
if(chunks==-1)	//using bp
{
if(bp[mark-1]<gbp1[j]-gene_buffer||chr[mark-1]<gchr[j]){break;}
}
else	//using weights
{
if(cumsum[mark-1]<=gbp1[j]||chr[mark-1]<gchr[j]){break;}
}
mark--;
}

while(mark<num_preds_use-1)	//see if we should move right
{
if(chunks==-1)	//using bp
{
if(bp[mark]>=gbp1[j]-gene_buffer||chr[mark+1]>gchr[j]){break;}
}
else	//using weights
{
if(cumsum[mark]>gbp1[j]||chr[mark+1]>gchr[j]){break;}
}
mark++;
}

//see if we are within - then find mark2 and set start and end for gene
if(chunks==-1)	//using bp
{
if(bp[mark]>=gbp1[j]-gene_buffer&&bp[mark]<=gbp2[j]+gene_buffer&&chr[mark]==gchr[j])
{
mark2=mark;
while(bp[mark2]<=gbp2[j]+gene_buffer&&chr[mark2]==gchr[j])
{
mark2++;
if(mark2==num_preds_use){break;}
}
gstarts[j]=mark;
gends[j]=mark2;
}
}
else	//using weights
{
if(cumsum[mark]>gbp1[j]&&cumsum[mark]<=gbp2[j]&&chr[mark]==gchr[j])
{
mark2=mark;
while(cumsum[mark2]<=gbp2[j]&&chr[mark2]==gchr[j])
{
mark2++;
if(mark2==num_preds_use){break;}
}
gstarts[j]=mark;
gends[j]=mark2;
}
}

if(gstarts[j]==-1&&chunks==-1&&chunksbp==-1)	//on chr but didn't find gene
{
if(misscount<10)
{printf("Warning, could not find any predictors within Gene %d (%s - Chr %d; %.0f - %.0f)\n", j+1, genenames[j], gchr[j], gbp1[j], gbp2[j]);}
misscount++;
}

}	//end of at correct chr
}	//end of loop through genes

if(misscount>0)
{printf("\nFor the chromosomes containing SNPs, %d genes contained no predictors\n\n", misscount);}


if(overlap==1)	//then when overlap, pick one
{
j=0;
while(j<num_genes_approx)
{
if(gstarts[j]!=-1){prev=j;break;}
j++;
}

for(j=prev+1;j<num_genes_approx;j++)
{
if(gstarts[j]!=-1)
{
if(gstarts[j]<gends[prev])	//then have an overlap (therefore must also be same chr)
{
mark=gstarts[j];
while(mark<gends[prev])
{
//see whether mark is closer to gene j, then gene prev, in which case goes in j, and break
if(chunks==-1)	//using bp
{
if(gbp1[j]-bp[mark]<bp[mark]-gbp2[prev]){break;}
}
else	//using weights
{
if(gbp1[j]-cumsum[mark]<cumsum[mark]-gbp2[prev]){break;}
}
mark++;
}
gends[prev]=mark;gstarts[j]=mark;
if(gstarts[prev]==gends[prev]){gstarts[prev]=-1;gends[prev]=-1;}
if(gstarts[j]==gends[j]){gstarts[j]=-1;gends[j]=-1;}
}
prev=j;
}
}
}	//end of overlap


//open file which stores gene details and write header lines
sprintf(filename1,"%sgene_details.txt",folder);
if((output1=fopen(filename1,"w"))==NULL)
{printf("Error writing to file %s - check write permission granted\n\n",filename1);exit(1);}

fprintf(output1, "Datafiles: %s %s\n", datafile, mapfile);

if(strcmp(kpredfile,"blank")==0){fprintf(output1, "Using All Predictors\n");}
else{fprintf(output1, "Using keepfile %s\n", kpredfile);}

fprintf(output1,"Gene_Number Partition Predictor_Start Predictor_End Gene_Name Gene_Length Gene_Weight Gene_Chr Gene_Start Gene_End Start_Predictor_BP End_Predictor_BP Min_Pvalue Min_Pvalue_Gene_BF\n");

count=0;count2=0;count3=0;
last=0;prev=0;prev2=0;
mark=-1;mark2=-1;
for(j=0;j<num_genes_approx;j++)
{
if(gstarts[j]!=-1&&(gstarts[j]!=mark||gends[j]!=mark2||(chunks!=-1||chunksbp!=-1)))	//skip if not found or boundaries match previous (unless considering genes)
{
mark=gstarts[j];mark2=gends[j];

//first get weightsum and min pvalue
weightsum=0;minpvalue=1;
for(j2=gstarts[j];j2<gends[j];j2++)
{
weightsum+=weights[j2];
if(pvalues[j2]<minpvalue){minpvalue=pvalues[j2];}
}
minpvalue2=1-pow(1-minpvalue,gends[j]-gstarts[j]);
if(minpvalue<1e-10){minpvalue2=minpvalue*(gends[j]-gstarts[j]);}

if(weightsum>minweight)
{
//see which predictors used, and add on number considered so far
for(j2=gstarts[j];j2<gends[j];j2++){usedpreds[j2]++;}
if(gstarts[j]<prev2){count2+=gends[j]-prev2;}
else{count2+=gends[j]-gstarts[j];}

//see which partition we should use (but don't allow gaps)
partition=(int)((count2-1)/partition_length)+1;
if(partition>last+1){partition=last+1;}
last=partition;

if(chunks!=-1||chunksbp!=-1)	//update genename
{sprintf(genenames[j],"Chunk_%d", count+1);}

if(chunks==-1)
{fprintf(output1, "%d %d %d %d %s %d %.2f %d %.0f %.0f %.0f %.0f ", count+1, partition, gstarts[j]+1, gends[j], genenames[j], gends[j]-gstarts[j], weightsum, gchr[j], gbp1[j], gbp2[j], bp[gstarts[j]], bp[gends[j]-1]);}
else
{fprintf(output1, "%d %d %d %d %s %d %f %d 0 0 %.0f %.0f ", count+1, partition, gstarts[j]+1, gends[j], genenames[j], gends[j]-gstarts[j], weightsum, gchr[j], bp[gstarts[j]], bp[gends[j]-1]);}

if(strcmp(pvafile,"blank")==0)
{fprintf(output1, "NA NA\n");}
else
{fprintf(output1, "%.2e %.2e\n", minpvalue, minpvalue2);}
count++;

if(gends[j]-gstarts[j]>count3){count3=gends[j]-gstarts[j];}

if((gstarts[j]<prev||gends[j]<prev2)&&bin==0)
{printf("Error, Gene %d (Predictors %d to %d) and Gene %d (Predictors %d to %d) are not monotonic, so data must be in a binary format (arguments \"--bfile\" or \"--speed\")\n\n", j, prev+1, prev2, j+1, gstarts[j]+1, gends[j]);exit(1);}

prev=gstarts[j];prev2=gends[j];
}
}
}	//end of j loop

fclose(output1);


//store the genic predictors and tallies
sprintf(filename2,"%sgene_preds.txt",folder);
sprintf(filename3,"%sgene_tallies.txt",folder);

if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to file %s - check write permission granted\n\n",filename2);exit(1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to file %s - check write permission granted\n\n",filename3);exit(1);}

count2=0;
for(j=0;j<num_preds_use;j++)
{
if(usedpreds[j]>0)
{fprintf(output2,"%d\n", keeppreds_use[j]+1);count2++;}
fprintf(output3,"%d\n", usedpreds[j]);
}
fclose(output2);fclose(output3);


printf("%d genes (or chunks) were found, spanning %d predictors, divided into %d partitions.\nThe details for those found are saved in %s\nThe longest gene contained %d predictors\n\n", count, count2, partition, filename1, count3);


free(cumsum);free(pvalues);free(usedpreds);

return(0);
}	//end of set_up_genes2

///////////////////////////

int get_genes_boundaries(char *folder, int partition, char **genenames, int *gchr, int *gstarts, int *gends, int *gpartitions, int num_genes, int num_preds_use)
{
int j;
int gmax;

int readint;
char readchar;
float readfloat;

char filename[500];
FILE *input;


//first open the details file - have alread checked present
sprintf(filename,"%sgene_details.txt", folder);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening file %s which should contain details of genes\nThis will have been created using argument \"--cut-genes\"\n\n", filename);exit(1);}

//skip the first three rows
readchar=0;
while(readchar!=10){(void)fscanf(input,"%c", &readchar);}
readchar=0;
while(readchar!=10){(void)fscanf(input,"%c", &readchar);}
readchar=0;
while(readchar!=10){(void)fscanf(input,"%c", &readchar);}

//now read in lines - index, partition, start, end, name
gmax=0;
for(j=0;j<num_genes;j++)
{
if(fscanf(input, "%d %d %d %d %s %d %f %d", &readint, gpartitions+j, gstarts+j, gends+j, genenames[j], &readint, &readfloat, gchr+j)!=8)
{printf("Error reading details for Gene %d (row %d of %s)\n\n", j+1, 3+j+1, filename);exit(1);}
gstarts[j]--;

//check sensible

if(gstarts[j]<0)
{printf("Error reading %s; the start predictor (%d) for Gene %d is less than 1; Has file been altered since creation using argument \"--cut-genes\"\n\n", filename, gstarts[j]+1, j+1);exit(1);}

if(gends[j]>num_preds_use&&num_preds_use!=-1)
{printf("Error reading %s; the end predictor (%d) for Gene %d is greater than total number of predictors (%d); Has file been altered since creation using argument \"--cut-genes\"\n\n", filename, gends[j]+1, j+1, num_preds_use);exit(1);}


//skip to end of row
readchar=0;
while(readchar!=10){(void)fscanf(input,"%c", &readchar);}

if(gpartitions[j]==partition&&gends[j]-gstarts[j]>gmax){gmax=gends[j]-gstarts[j];}
}	//end of j loop

fclose(input);

if(gmax==0&&partition!=-1)
{printf("Error reading %s; there are no genes (of non-zero length) in Partition %d\n\n", filename, partition);exit(1);}

return(gmax);
}	//end of get_genes_boundaries

///////////////////////////

int join_reml(char *folder, char **genenames, int *gstarts, int *gends, int *gpartitions, int *gchr, int num_genes, char **prednames, float cut1, float cut2, int scoretest)
{
int j, g, partition, count, count2;
int num_parts, mpheno, start, end, new1, new2;

int readint, readphen;
char readchar, readstring[100];

int *usedgenes;
float values[11], *gdetails, *pvalues;

FILE *input1, *output2, *output3, *output4;
char filename1[500], filename2[500], filename3[500], filename4[500];


gdetails=malloc(sizeof(float)*num_genes*11);
usedgenes=malloc(sizeof(int)*num_genes);
for(g=0;g<num_genes;g++){usedgenes[g]=0;}

//find how many partitions
num_parts=0;
for(g=0;g<num_genes;g++)
{
if(gpartitions[g]>num_parts){num_parts=gpartitions[g];}
}

//open each file in turn - keep track of which phen was tested
mpheno=0;
for(partition=0;partition<num_parts;partition++)
{
//check partition results file present and get size
sprintf(filename1, "%sregress%d", folder, partition+1);
if((input1=fopen(filename1,"r"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
fclose(input1);
count=countrows(filename1)-1;

//skip header row
if((input1=fopen(filename1,"r"))==NULL)
{printf("Error opening %s\n\n", filename1);exit(1);}
if(fscanf(input1, "%s ", readstring)!=1)
{printf("Error reading first element of %s\n\n", filename1);exit(1);}
if(strcmp(readstring,"Gene_Number")!=0)
{printf("Error reading %s; first element should be \"Gene_Number\" (not %s)\n\n", filename1, readstring);exit(1);}
readchar=0;
while(readchar!=10){(void)fscanf(input1,"%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input1, "%d %s %d %f %f %f %f %f %f %f %f %f %f %f ", &readint, readstring, &readphen, values, values+1, values+2, values+3, values+4, values+5, values+6, values+7, values+8, values+9, values+10)!=14)
{printf("Error reading results from Row %d of %s\n\n", j+2, filename1);exit(1);}

if(readint>num_genes||readint<1)
{printf("Error reading %s; file not consistent with details file created by \"--cut-genes\"\n\n", filename1);exit(1);}
readint--;
usedgenes[readint]++;

if(strcmp(genenames[readint],readstring)!=0)
{printf("Error Doug\n\n");exit(1);}

if(mpheno==0){readphen=mpheno;}
if(readphen!=mpheno)
{printf("Error reading %s; contains results for Phenotype %d, but earlier results were for Phenotype %d\n\n", filename1, readphen, mpheno);exit(1);}

gdetails[readint]=values[0];gdetails[readint+num_genes]=values[1];gdetails[readint+2*num_genes]=values[2];
gdetails[readint+3*num_genes]=values[3];gdetails[readint+4*num_genes]=values[4];
gdetails[readint+5*num_genes]=values[5];gdetails[readint+6*num_genes]=values[6];
gdetails[readint+7*num_genes]=values[7];gdetails[readint+8*num_genes]=values[8];
gdetails[readint+9*num_genes]=values[9];gdetails[readint+10*num_genes]=values[10];
}	//end of j loop
fclose(input1);
}	//end of partition loop


//check all present
for(g=0;g<num_genes;g++)
{
if(usedgenes[g]==0){printf("Error, results not found for Gene %d (%s)\n\n", g+1, genenames[g]);exit(1);}
if(usedgenes[g]>1){printf("Error, multiple results found for Gene %d (%s)\n\n", g+1, genenames[g]);exit(1);}
}

//now store
sprintf(filename2, "%sregressALL", folder);
if((output2=fopen(filename2, "w"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
fprintf(output2, "Gene_Number Gene_name Phen_Number REML_Her REML_SD LRT_Stat LRT_P Score_Stat Score_P MA_Delta MA_SD REML_BF Gene_Length Gene_Weight\n");
for(g=0;g<num_genes;g++)
{fprintf(output2, "%d %s %d %.6f %.6f %.2f %.4e %.2f %.4e %.4f %.4f %.2f %d %.2f\n", g+1, genenames[g], mpheno, gdetails[g], gdetails[g+num_genes], gdetails[g+2*num_genes], gdetails[g+3*num_genes], gdetails[g+4*num_genes], gdetails[g+5*num_genes], gdetails[g+6*num_genes], gdetails[g+7*num_genes], gdetails[g+8*num_genes], (int)gdetails[g+9*num_genes], gdetails[g+10*num_genes]);}
fclose(output2);

if(cut1!=-1)	//prepare regions for MultiBLUP
{
pvalues=gdetails+3*num_genes;
if(scoretest==1){pvalues=gdetails+5*num_genes;}

count=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(pvalues[g]<=cut1)	//gene is sig
{
count++;
pvalues[g]=2;
start=g;end=g;

//move left
while(1)
{
new1=start-1;new2=start-2;
if(new1<0){break;}
if(new2<0){new2=0;}
if(gchr[new1]!=gchr[g]){break;}
if(gchr[new2]!=gchr[g]){new2=new1;}
if(gends[new2]<gstarts[new1]){new2=new1;}
if(gends[new1]<gstarts[start]){break;}

if(pvalues[new1]>cut2&&pvalues[new2]>cut2){break;}
pvalues[new1]=2;
start=new1;
}

//move right
while(1)
{
new1=end+1;new2=end+2;
if(new1>=num_genes){break;}
if(new2>=num_genes){new2=num_genes-1;}
if(gchr[new1]!=gchr[g]){break;}
if(gchr[new2]!=gchr[g]){new2=new1;}
if(gstarts[new2]>gends[new1]){new2=new1;}
if(gstarts[new1]>gends[end]){break;}

if(pvalues[new1]>cut2&&pvalues[new2]>cut2){break;}
pvalues[new1]=2;
end=new1;
}

sprintf(filename3, "%sregion%d", folder, count);
if((output3=fopen(filename3, "w"))==NULL)
{printf("Error opening %s\n\n", filename3);exit(1);}
for(j=gstarts[start];j<gends[end];j++)
{fprintf(output3,"%s\n",prednames[j]);count2++;}
fclose(output3);
}	//end of gene g sig
}	//end of g loop

printf("For MultiBLUP, found %d regions spanning %d predictors (sig thresholds %f and %f)\n", count, count2, cut1, cut2);
sprintf(filename4, "%sregion_number", folder);
if((output4=fopen(filename4, "w"))==NULL)
{printf("Error opening %s\n\n", filename4);exit(1);}
fprintf(output4, "%d\n", count);
fclose(output4);
}		//end of preparing regions

free(gdetails);free(usedgenes);

return(0);
}	//end of join_reml

///////////////////////////



