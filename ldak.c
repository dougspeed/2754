//??? - sorting strings once only?
//??? at some point make a guess of memory requirements?
//??? read eigenfile
//??? get bf for generalized reml
//??? need to look at role of ea in bivariate

/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//to compile this code:

//if you have the gnu linear programming (glpk) and scientific (gsl) libraries installed, use
//gcc source/ldak.c -lm -lglpk -lgsl -lblas -llapack -o ldak.out
//(note, you may have to toggle which of glpk/glpk.h and glpk.h is linked to below)

//if you do not have these installed, you can try including the files in the download folder
//gcc source/ldak.c -lc libraries/libgsl.a libraries/libglpk.a libraries/liblapack.a libraries/libblas.a libraries/libc.a libraries/libgfortran.a libraries/libm.a -I source/ldak -static -o ldak.out 


//library includes

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<errno.h>
#include<time.h>


//use one of these two includes if you have glpk installed on your computer
//#include<glpk.h>
//#include<glpk/glpk.h>

//and use these includes if you have gsl installed on your computer
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>


//will export functions from c

extern void dgemm_();
extern void dgemv_();
extern void dpotrf_();
extern void dpotri_();
//extern void dpotrs_();
//extern void dgetrf_();
//extern void dgetrs_();
extern void dsyev_();


#include "sort.c"
#include "permute.c"
#include "norm.c"
#include "fileops.c"
#include "weights.c"
#include "kinships.c"
#include "readdata.c"
#include "genes.c"
#include "remlgenes.c"
#include "remlbivar.c"
#include "remlmulti.c"
#include "blup.c"


int main (int argc, const char * argv[])
{
//set random variable seed
srand((unsigned)time(NULL));
//set seed for normal generator
zigset(rand());

//and something for gsl random num distributors
//const gsl_rng_type * gslT;
//gsl_rng * gslr;
//gsl_rng_env_setup();
//gslT = gsl_rng_default;
//gslr = gsl_rng_alloc (gslT);

printf("\nLDAK - Software for obtaining Linkage Disequilibrium Adjusted Kinship estimates and Loads More\n");
printf("Help pages at http://dougspeed.com/ldak\n\n");

//declare variables
#include "declare.c"

//deal with command line arguments (come in pairs)
#include "readargs.c"

//manipulate some file names
#include "append.c"

//do some checks that files are present when provided or required
#include "checks1.c"

//and some more mode-specific checks
#include "checks2.c"

//and now check the details files for weights, kinships or genes
if(mode==2||mode==3||mode==5||mode==6||mode==8||mode==9||mode==10||mode==11||mode==12)
{
#include "checks3.c"
}

/////////

//see if changing famfile (only allowed if not reading data)

if((mode==16&&num_regs==0&&num_kins>0)||mode==17||mode==18)	//get from kinstems[0]
{printf("Will set ids based on %s.grm.id\n", kinstems[0]);
sprintf(famfile,"%s.grm.id", kinstems[0]);}

if(mode==16&&num_kins==0&&num_regs==0)	//get from respfile
{printf("Will set ids based on %s\n", respfile);
sprintf(famfile,"%s", respfile);}

//sort which samples and predictors we are using - also deals with re-calc weights and subtracting regions
#include "getnums.c"

if(num_samples!=-1)
{
printf("Original number of samples: %d --- Number being used: %d\n", num_samples, num_samples_use);
if(num_samples_use<3){printf("Can't continue with so few samples\n\n");exit(1);}
}
if(num_preds!=-1)
{printf("Original number of predictors: %d --- Number being used: %d\n", num_preds, num_preds_use);}
printf("\n");

///////////////////////////
///////////////////////////

//print out parameters depending on mode
#include "param.c"

//set data_length to num_preds_use and keeppreds_use to keeppreds, then change when necessary
//data_length always refers to main predictors (not regions)
if(num_preds!=-1)
{
data_length=num_preds_use;
keeppreds_use=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){keeppreds_use[j]=keeppreds[j];}
}

/////////

if(mode==8||mode==9||mode==10||mode==11||mode==12)	//will require gene details
{
sprintf(filename, "%sgene_details.txt", folder);
num_genes=countrows(filename)-3;

genenames=malloc(sizeof(char*)*num_genes);
for(j=0;j<num_genes;j++){genenames[j]=malloc(sizeof(char)*100);}
gchr=malloc(sizeof(int)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
gpartitions=malloc(sizeof(int)*num_genes);

gmax=get_genes_boundaries(folder, partition, genenames, gchr, gstarts, gends, gpartitions, num_genes, num_preds_use);
}

if(num_regs>0&&(mode==9||mode==16||mode==24))	//will be using regions
{
rkeeppreds_use=malloc(sizeof(int)*num_preds_use);
regindex=malloc(sizeof(int*)*num_regs);
rdata_length=read_regions(regpref, num_regs, mapfile, regindex, rkeeppreds_use);
printf("In total, there are %d (unique) region predictors\n\n", rdata_length);
}

/////////

//if(mode==1)	//cut for weights - default fine
if(mode==2)	//calc weights, get boundaries
{
if(strcmp(weightsfile,"blank")==0)
{data_length=get_section_boundaries(folder, section, keeppreds_use, keeppreds, num_preds_use, 0);}
else
{data_length=get_section_boundaries(folder, section, keeppreds_use, keeppreds, num_preds_use, 1);}
}
//if(mode==3)	//join weights - default fine

/////////

//if(mode==4)	//cut for kins - default fine
if(mode==5)	//calc kins, get boundaries
{
data_length=get_partition_boundaries(folder, partition, keeppreds_use, keeppreds, num_preds_use, mapfile);
}
//if(mode==6)	//join kins - default fine

/////////

//if(mode==7)	//cut for genes - default fine
//if(mode==8)	//gene kins - default + gene details fine
//if(mode==9)	//gene reml - default + gene details fine
//if(mode==10)	//join gene reml - default fine
//if(mode==11)	//gene bivar - default + gene details fine
//if(mode==12)	//join gene bivar - default fine

/////////

//if(mode==16)	//generalized reml - default fine
//if(mode==17||mode==18)	//decompose kinships - default fine

if(mode==19||mode==20)	//calc blups or scores - get predictor effects
{
allpreds=malloc(sizeof(char*)*num_preds_use);
for(j=0;j<num_preds_use;j++){allpreds[j]=malloc(sizeof(char)*100);}
alla1=malloc(sizeof(char)*num_preds_use);
alla2=malloc(sizeof(char)*num_preds_use);
read_mapfile(mapfile, NULL, allpreds, NULL, alla1, alla2, num_preds_use, keeppreds);

allcentres=malloc(sizeof(float*)*(num_kins+num_regs));
allfactors=malloc(sizeof(float*)*(num_kins+num_regs));
wsums=malloc(sizeof(float)*num_kins);

data_length=set_up_blup1(allcentres, allfactors, wsums, num_kins, kinstems, num_regs, regfile, allpreds, alla1, alla2, num_preds_use, keeppreds, keeppreds_use, famfile, adjust);

//for regions, values in allfactors are effects
effects=malloc(sizeof(double*)*(num_kins+num_regs+1));
for(k=0;k<num_kins+num_regs+1;k++){effects[k]=malloc(sizeof(double)*data_length);}
for(r=0;r<num_regs;r++)
{
for(j=0;j<data_length;j++){effects[num_kins+r][j]=allfactors[num_kins+r][j];}
}
for(j=0;j<data_length;j++){effects[num_kins+num_regs][j]=0;}
for(j=0;j<num_preds_use;j++){free(allpreds[j]);}free(allpreds);
free(alla1);free(alla2);
}	//end of prepare for blup

////////

if(mode==21)	//make phens ???
{}
if(mode==22)	//make snps ???
{}

//if(mode==23)	//add kins - default fine
//if(mode==24&&num_regs==0)	//subtract kins - default fine
if(mode==24&&num_regs>0)	//get (max) number of predictors
{
sprintf(filename,"%s.grm.details",kinstems[0]);
num_predsb=countrows(filename)-1;
}

//if(mode==25||mode==26||mode==27)	//converting data types - default fine

///////////////////////////
///////////////////////////

if(num_samples!=-1)	//get ids
{
ids1=malloc(sizeof(char*)*num_samples_use);ids2=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++){ids1[i]=malloc(sizeof(char)*100);ids2[i]=malloc(sizeof(char)*100);}
read_famfile(famfile, ids1, ids2, num_samples_use, keepsamps);
}

if(num_preds!=-1)	//read in regular (or gene) predictor details (unnecessary for some modes)
{
chr=malloc(sizeof(int)*data_length);
prednames=malloc(sizeof(char*)*data_length);
for(j=0;j<data_length;j++){prednames[j]=malloc(sizeof(char)*100);}
bp=malloc(sizeof(double)*data_length);
al1=malloc(sizeof(char)*data_length);
al2=malloc(sizeof(char)*data_length);

if(strcmp(chiamofile,"blank")!=0)
{read_mapfile(mapfile, chr, prednames, bp, NULL, NULL, data_length, keeppreds_use);}
if(strcmp(chiamofile,"blank")==0)
{read_mapfile(mapfile, chr, prednames, bp, al1, al2, data_length, keeppreds_use);}
}

///////////////////////////

if(mode==5||mode==7||mode==8||mode==9||mode==11)	//sort out regular (or gene) predictor weights
{
weights=malloc(sizeof(float)*data_length);

if(ignoreweights==1)	//set weights to 1
{
for(j=0;j<data_length;j++){weights[j]=1;}
}
else	//read weights from weightsfile
{
read_weightfile(weightsfile, weights, data_length, keeppreds_use, prednames, mapfile);
printf("First few weights are:\n");
for(j=0;j<data_length;j++){if(j<5){printf("%f   ", weights[j]);}}
printf("\n\n");
}
}

/////////

if(mode==9||mode==11||mode==16||(mode==17&&num_resps>0))	//read resp and covar (might update ids)
{
resp=malloc(sizeof(double)*num_samples_use*num_resps_use);
respindex=malloc(sizeof(int*)*num_resps_use);
for(m=0;m<num_resps_use;m++){respindex[m]=malloc(sizeof(int)*(num_samples_use+1));}

read_respfile(respfile, resp, ids1, ids2, num_samples_use, num_resps_use, keepresps, num_resps, missingvalue, respindex);

//update details based on which ids will be used
num_samples_use=update_resps(resp, respindex, ids1, ids2, keepsamps, num_samples_use, num_resps_use);

if(permute==1)	//permute
{
order=malloc(sizeof(int)*num_samples_use);
resptemp=malloc(sizeof(double)*num_samples_use);
for(i=0;i<num_samples_use;i++){order[i]=i;}
for(m=0;m<num_resps_use;m++)
{
permute_int(order, num_samples_use);
for(i=0;i<num_samples_use;i++){resptemp[i]=resp[i+m*num_samples_use];}
for(i=0;i<num_samples_use;i++){resp[i]=resptemp[order[i]+m*num_samples_use];}
}
free(order);free(resptemp);
}

covar=malloc(sizeof(double)*num_samples_use*num_covars);
for(i=0;i<num_samples_use;i++){covar[i]=1;}

if(strcmp(covarfile,"blank")!=0)
{read_covarfile(covarfile, covar+num_samples_use, ids1, ids2, num_samples_use, num_covars-1, missingvalue);}

//fill Y and Z for response 1
Ya=malloc(sizeof(double)*num_samples_use);
Za=malloc(sizeof(double)*num_samples_use*num_covars);
for(i=0;i<respindex[0][0];i++)
{
Ya[i]=resp[respindex[0][1+i]];
for(j=0;j<num_covars;j++){Za[i+j*respindex[0][0]]=covar[respindex[0][1+i]+j*num_samples_use];}
}
if(mode==11)	//and for response 2
{
Yb=malloc(sizeof(double)*num_samples_use);
Zb=malloc(sizeof(double)*num_samples_use*num_covars);
for(i=0;i<respindex[1][0];i++)
{
Yb[i]=resp[respindex[1][1+i]];
for(j=0;j<num_covars;j++){Zb[i+j*respindex[1][0]]=covar[respindex[1][1+i]+j*num_samples_use];}
}
}
}	//end of reading resp and covar

/////////

if(num_kins>0&&(mode==9||mode==16||mode==17))	//read kins
{
mkins=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++){mkins[k]=malloc(sizeof(double*)*num_samples_use*num_samples_use);}

for(k=0;k<num_kins;k++)
{read_kinfile(kinstems[k], mkins[k], NULL, num_samples_use, ids1, ids2, 0);}
}

///////////////////////////

if(num_regs>0&&(mode==9||mode==16||mode==24))	//get region predictors and details
{
rprednames=malloc(sizeof(char*)*rdata_length);
for(j=0;j<rdata_length;j++){rprednames[j]=malloc(sizeof(char)*100);}
ral1=malloc(sizeof(char)*rdata_length);
ral2=malloc(sizeof(char)*rdata_length);

if(strcmp(chiamofile,"blank")!=0)
{read_mapfile(mapfile, NULL, rprednames, NULL, NULL, NULL, rdata_length, rkeeppreds_use);}
if(strcmp(chiamofile,"blank")==0)
{read_mapfile(mapfile, NULL, rprednames, NULL, ral1, ral2, rdata_length, rkeeppreds_use);}

rweights=malloc(sizeof(float)*rdata_length);

if(ignoreweights==1)	//set weights to 1
{
for(j=0;j<rdata_length;j++){rweights[j]=1;}
}
else	//read weights from weightsfile
{
read_weightfile(weightsfile, rweights, rdata_length, rkeeppreds_use, rprednames, mapfile);
printf("First few region weights are:\n");
for(j=0;j<rdata_length;j++){if(j<5){printf("%f   ", rweights[j]);}}
printf("\n\n");
}

rdata=malloc(sizeof(double)*num_samples_use*rdata_length);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

rkeepcentres=malloc(sizeof(double)*rdata_length);
rkeepmults=malloc(sizeof(double)*rdata_length);

(void) read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, 0, rdata, num_samples_use, keepsamps, 0, rdata_length, rkeeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, ral1, ral2, -1);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}

//get stats but don't standardize
(void)get_stats_qc(rkeepcentres, rkeepmults, rdata, num_samples_use, rdata_length, minmaf, maxmaf, minvar, minobs, missingvalue, power, rprednames, 0, rweights);

if(mode==9||mode==16)	//prune and update stats
{prune_regions(num_regs, regindex, rdata, num_samples_use, rdata_length, prune, rkeepmults, rweights);}

//now standardize
(void)get_stats_qc(rkeepcentres, rkeepmults, rdata, num_samples_use, rdata_length, minmaf, maxmaf, minvar, minobs, missingvalue, power, rprednames, 3, rweights);

//get an approx max
rmax=0;for(r=0;r<num_regs;r++){rmax+=regindex[r][0];}
}

///////////////////////////

if(mode==9||mode==16)	//set shortcut	??? think about this a bit more
{
count=0;
for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++)
{
j2=regindex[r][1+j];
if(rkeepmults[j2]!=-1&&rweights[j2]>0){count++;}
}
}
if(count<=num_samples_use&&num_kins<2){shortcut=1;}
if(count>num_samples_use&&strcmp(eigenfile,"blank")!=0)
{printf("Warning, there are more region predictors (%d) than samples (%d), so will not be using eigen-decomposition\n\n", count, num_samples_use);}

//shortcut=0;
//printf("Shortcut is %d\n\n", shortcut);
}

////////

if(strcmp(eigenfile,"blank")!=0&&shortcut==1)	//read decompostion
{
printf("Reading eigen-decomposition from %s\n", eigenfile);
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);
read_eigenfile(eigenfile, U, E, num_samples_use, ids1, ids2);
}

if(mode==17||mode==18||(shortcut==1&&num_kins==1&&strcmp(eigenfile,"blank")==0))	//perform decomposition
{
printf("Performing eigen-decomposition\n");
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);

for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<num_samples_use;i2++){U[i+i2*num_samples_use]=mkins[0][i+i2*num_samples_use];}
}
lwork=-1;
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, &wkopt, &lwork, &info );
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, work, &lwork, &info);
free(work);
printf("Eigen-decomposition complete\n\n");
}

////////

if(mode==19&&num_kins>0)	//do remainder of blup prep (which required ids)
{
mG2=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++){mG2[k]=malloc(sizeof(double)*num_samples_use);}
set_up_blup2(mG2, blupfile, ids1, ids2, num_samples_use, num_kins, num_regs);
}

///////////////////////////
///////////////////////////

//perform work for mode

if(mode==1)	//cut for weights
{
if(strcmp(weightsfile,"blank")==0)
{set_up_sections(folder, chr, bp, num_preds_use, section_length, section_buffer, mapfile, datafile, bpredfile, 0);}
else
{set_up_sections(folder, chr, bp, num_preds_use, section_length, section_buffer, mapfile, datafile, bpredfile, 1);}
}

if(mode==2)	//calc cor then solve linear programming problem
{
//read in data
data=malloc(sizeof(double)*num_samples_use*data_length);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

keepcentres=malloc(sizeof(double)*data_length);
keepmults=malloc(sizeof(double)*data_length);

(void) read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, 0, data, num_samples_use, keepsamps, 0, data_length, keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}

(void)get_stats_qc(keepcentres, keepmults, data, num_samples_use, data_length, minmaf, maxmaf, minvar, minobs, missingvalue, power, prednames, 1, NULL);

//calc correlations

//solve simplex - if fails, suggest fudge

free(data);free(keepcentres);free(keepmults);
}

if(mode==3)	//join weights
{
}

/////////

if(mode==4)	//cut for kins
{
set_up_partitions(folder, num_preds_use, partition_length, partpref, num_parts, chr, bychr, prednames, datafile, mapfile, bpredfile);
}

if(mode==5)	//calc kinship
{
//will read in data a bit at a time
data=malloc(sizeof(double)*num_samples_use*bitsize);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

keepcentres=malloc(sizeof(double)*data_length);
keepmults=malloc(sizeof(double)*data_length);

kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
for(i=0;i<num_samples_use*num_samples_use;i++){kins[i]=0;}

bittotal=(int)((data_length-1)/bitsize+1);
current=0;weightsum=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating kinships for Chunk %d of %d in Partition %d\n", bit+1, bittotal, partition);

current=read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, current, data, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);

weightsum+=get_stats_qc(keepcentres+bitstart, keepmults+bitstart, data, num_samples_use, bitlength, minmaf, maxmaf, minvar, minobs, missingvalue, power, prednames, 2, weights+bitstart);

alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);
}

if(weightsum==0)
{printf("Error, no predictors remain in Partition %d\n\n", partition);exit(1);}

//save kins
sprintf(outfile, "%skinship%d", folder, partition);
write_kins(outfile, kins, NULL, weightsum, ids1, ids2, num_samples_use, keeppreds_use, prednames, keepcentres, keepmults, weights, al1, al2, data_length, kinraw, kingz);

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}
free(data);free(keepcentres);free(keepmults);free(kins);
}

if(mode==6)	//join kinships
{
manip_kins(outfile, num_kins, kinstems, samedata, kindetails, kinraw, kingz, 0);
}

/////////

if(mode==7)	//cut genes
{
num_genes_approx=set_up_genes1(genefile, chunksbp, chunks, chr, bp, weights, num_preds_use)+10;

genenames=malloc(sizeof(char*)*num_genes_approx);
for(j=0;j<num_genes_approx;j++){genenames[j]=malloc(sizeof(char)*100);}
gchr=malloc(sizeof(int)*num_genes_approx);
gbp1=malloc(sizeof(double)*num_genes_approx);
gbp2=malloc(sizeof(double)*num_genes_approx);
gstarts=malloc(sizeof(int)*num_genes_approx);
gends=malloc(sizeof(int)*num_genes_approx);
for(j=0;j<num_genes_approx;j++){gstarts[j]=-1;gends[j]=-1;}

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{set_up_genes2(folder, gchr, gbp1, gbp2, genenames, gstarts, gends, genefile, chunksbp, chunks, chr, bp, weights, num_preds_use, gene_buffer, minweight, overlap, partition_length, pvafile, keeppreds_use, mapfile, datafile, bpredfile, 0);}
else
{set_up_genes2(folder, gchr, gbp1, gbp2, genenames, gstarts, gends, genefile, chunksbp, chunks, chr, bp, weights, num_preds_use, gene_buffer, minweight, overlap, partition_length, pvafile, keeppreds_use, mapfile, datafile, bpredfile, 1);}

for(j=0;j<num_genes_approx;j++){free(genenames[j]);}free(genenames);
free(gchr);free(gbp1);free(gbp2);free(gstarts);free(gends);
}

/////////

if(mode==8)	//get gene-based kinships
{
//will read in data a gene at a time
data=malloc(sizeof(double)*num_samples_use*gmax);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

keepcentres=malloc(sizeof(double)*data_length);
keepmults=malloc(sizeof(double)*data_length);

kins=malloc(sizeof(double)*num_samples_use*num_samples_use);

//count out how many genes we will consider
total=0;for(g=0;g<num_genes;g++){if(gpartitions[g]==partition){total++;}}

current=0;prev1=0;prev2=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(gpartitions[g]==partition)
{
printf("Calculating kinships for Gene %d out of %d in Partition %d\n", count2+1, total, partition);

glength=gends[g]-gstarts[g];
count=0;weightsum=0;

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)	//shuffle along data
{
for(j=gstarts[g];j<prev2;j++)
{
for(i=0;i<num_samples_use;i++)
{data[i+count*num_samples_use]=data[i+(j-prev1)*num_samples_use];}
weightsum+=weights[gstarts[g]+count];
count++;
}
}

current=read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, current, data+count*num_samples_use, num_samples_use, keepsamps, gstarts[g]+count, gends[g], keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);

weightsum+=get_stats_qc(keepcentres+gstarts[g]+count, keepmults+gstarts[g]+count, data+count*num_samples_use, num_samples_use, gends[g]-gstarts[g]-count, minmaf, maxmaf, minvar, minobs, missingvalue, power, prednames, 2, weights+gstarts[g]+count);

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &glength, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);

sprintf(outfile, "%sgeneship%d", folder, g+1);
write_kins(outfile, kins, NULL, weightsum, ids1, ids2, num_samples_use, keeppreds_use+gstarts[g], prednames+gstarts[g], keepcentres+gstarts[g], keepmults+gstarts[g], weights+gstarts[g], al1+gstarts[g], al2+gstarts[g], glength, kinraw, kingz);

prev1=gstarts[g];prev2=gends[g];count2++;
}	//end of considering gene g
}	//end of g loop

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}
free(data);free(keepcentres);free(keepmults);free(kins);
}	//end of mode=8

/////////

if(mode==9||mode==11)	//gene-based testing
{
//will read in data a gene at a time
data=malloc(sizeof(double)*num_samples_use*gmax);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

keepcentres=malloc(sizeof(double)*data_length);
keepmults=malloc(sizeof(double)*data_length);

Xstarts=malloc(sizeof(int)*(num_regs+1));Xends=malloc(sizeof(int)*(num_regs+1));
Xrec=malloc(sizeof(int)*(rmax+gmax));Xrev=malloc(sizeof(int)*(rmax+gmax));
Xsums=malloc(sizeof(double)*(num_regs+1));

//prep for first response - might have regions
Xa=malloc(sizeof(double)*respindex[0][0]*(rmax+gmax));
if(num_regs+num_kins==0)
{
ZTYa=malloc(sizeof(double)*num_covars);
ZTZa=malloc(sizeof(double)*num_covars*num_covars);
prepare_gene_reml(Ya, Za, respindex[0][0], num_covars, &YTYa, ZTYa, ZTZa, &detZTZa, &YTCYa);
}
else
{
Mnull=malloc(sizeof(double)*respindex[0][0]*respindex[0][0]);
fill_X(Xa, Xstarts, Xends, Xrec, Xrev, Xsums, respindex[0], num_samples_use, num_regs, rdata, regindex, rkeepmults, rweights, data, 0, NULL, NULL, &adjust2);
multi_reml(Ya, Za, Xa, respindex[0][0], num_covars, num_regs, 0, Xstarts, Xends, Xrec, Xrev, Xsums, num_kins, mkins, shortcut, U, E, remla, NULL, NULL, NULL, rdata_length, rkeepcentres, rkeepmults, rweights, rprednames, ral1, ral2, NULL, NULL, NULL, 0, 0, 0, Mnull);
likenull=remla[2];
}

if(mode==11)	//prep for second response (at moment will always have num_kins+num_regs=0)
{
Xb=malloc(sizeof(double)*respindex[1][0]*(rmax+gmax));
ZTYb=malloc(sizeof(double)*num_covars);
ZTZb=malloc(sizeof(double)*num_covars*num_covars);
prepare_gene_reml(Yb, Zb, respindex[1][0], num_covars, &YTYb, ZTYb, ZTZb, &detZTZb, &YTCYb);
}

//open results files
if(mode==9)
{
sprintf(filename,"%sregress%d", folder, partition);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s\n\n",filename);exit(1);}
fprintf(output, "Gene_Number Gene_Name Phen_Number REML_Her REML_SD Likelihood LRT_Stat LRT_P Score_Stat Score_P MA_Delta MA_SD REML_BF Gene_Length Gene_Weight\n");

sprintf(filename2,"%seffects%d", folder, partition);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s\n\n", filename2);exit(1);}
fprintf(output2, "Gene_Number Predictor Effect\n");
fclose(output2);
}
if(mode==11)
{
sprintf(filename2,"%seffects%d", folder, partition);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s\n\n", filename2);exit(1);}
fprintf(output2, "Gene_Number Predictor Effect\n");

sprintf(filename,"%sregressA%d", folder, partition);
if((output2=fopen(filename,"w"))==NULL)
{printf("Error writing to %s\n\n", filename);exit(1);}
fprintf(output, "Gene_Number Gene_Name Phen_Number REML_HerA REML_SDA LikelihoodA LRT_StatA LRT_PA Score_StatA Score_PA MA_DeltaA MA_SDA REML_BFA Gene_Length Gene_Weight\n");

sprintf(filename2,"%sregressB%d", folder, partition);
if((output3=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s\n\n", filename2);exit(1);}
fprintf(output2, "Gene_Number Gene_Name Phen_Number REML_HerB REML_SDB LikelihoodB LRT_StatB LRT_P Score_StatB Score_PB MA_DeltaB MA_SDB REML_BFB Gene_Length Gene_Weight\n");

sprintf(filename3,"%sbivar%d", folder, partition);
if((output=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s\n\n", filename3);exit(1);}
fprintf(output3, "Gene_Number Gene_Name Correlation Cor_SD REML_Her1 REML_Her2 LRT_Stat_0 LRT_P_0 Score_Stat_0 Score_Pva_0 Gene_Length Gene_Weight\n");
}

//count out how many genes we will consider
total=0;for(g=0;g<num_genes;g++){if(gpartitions[g]==partition){total++;}}

current=0;prev1=0;prev2=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(gpartitions[g]==partition)
{
if(count2%100==0||(count2%10==0&&num_kins+num_regs>0))
{
printf("Testing Gene %d out of %d in Partition %d\n", count2, total, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL){printf("Error re-opening %s\n\n",filename);exit(1);}
if(mode==11)
{
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL){printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL){printf("Error re-opening %s\n\n",filename3);exit(1);}
}
}

glength=gends[g]-gstarts[g];
count=0;
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)	//shuffle along data
{
for(j=gstarts[g];j<prev2;j++)
{
for(i=0;i<num_samples_use;i++)
{data[i+count*num_samples_use]=data[i+(j-prev1)*num_samples_use];}
count++;
}
}

current=read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, current, data+count*num_samples_use, num_samples_use, keepsamps, gstarts[g]+count, gends[g], keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);

get_stats_qc(keepcentres+gstarts[g]+count, keepmults+gstarts[g]+count, data+count*num_samples_use, num_samples_use, glength-count, minmaf, maxmaf, minvar, minobs, missingvalue, power, prednames, 2, weights+gstarts[g]+count);

//calc gene-based reml for first response
glength2=fill_X(Xa, Xstarts, Xends, Xrec, Xrev, Xsums, respindex[0], num_samples_use, num_regs, rdata, regindex, rkeepmults, rweights, data, glength, keepmults, weights, &adjust2);

if(num_kins+num_regs==0)	//easy case
{gene_reml(remla, Ya, Za, Xa, respindex[0][0], num_covars, glength2, Xsums[0], YTYa, ZTYa, ZTZa, detZTZa, YTCYa, priora, priorb);}
else	//general case
{multi_reml(Ya, Za, Xa, respindex[0][0], num_covars, num_regs, g+1, Xstarts, Xends, Xrec, Xrev, Xsums, num_kins, mkins, shortcut, U, E, remla, filename2, NULL, NULL, rdata_length, rkeepcentres, rkeepmults, rweights, rprednames, ral1, ral2, keepmults+gstarts[g], weights+gstarts[g], prednames+gstarts[g], adjust, adjust2, likenull, Mnull);}

fprintf(output, "%d %s %d %.6f %.6f %.2f %.2f %.4e %.2f %.4e %.4f %.4f %.2f %d %.2f\n", g+1, genenames[g], keepresps[0]+1, remla[0], remla[1], remla[2], remla[3], remla[4], remla[5], remla[6], remla[7], remla[8], remla[9], glength, Xsums[num_regs]);

if(mode==11)	//calc gene-based bivariate (currently must have num_kins+num_regs=0)
{
if(fill_X(Xb, NULL, NULL, NULL, NULL, NULL, respindex[1], num_samples_use, num_regs, rdata, regindex, rkeepmults, rweights, data, glength, keepmults, weights, &adjust2)!=glength2)
{printf("Doug error 5\n");}

gene_reml(remlb, Yb, Zb, Xb, respindex[1][0], num_covars, glength2, Xsums[0], YTYb, ZTYb, ZTZb, detZTZb, YTCYb, priora, priorb);

fprintf(output2, "%d %s %d %.6f %.6f %.2f %.2f %.4e %.2f %.4e %.4f %.4f %.2f %d %.2f\n", g+1, genenames[g], keepresps[0]+1, remlb[0], remlb[1], remlb[2], remlb[3], remlb[4], remlb[5], remlb[6], remlb[7], remlb[8], remlb[9], glength, Xsums[num_regs]);

//set bivar - gene sizes, rho, SDrho, herA, SDA, herB, SDB, LRTstat, LRTPval, scorestat, scorePval, eA, sigA, eB, sigB
bivar[0]=0;bivar[1]=-1;bivar[2]=remla[0];bivar[3]=remla[1];bivar[4]=remlb[0];bivar[5]=remlb[1];
bivar[6]=0;bivar[7]=1;bivar[8]=0;bivar[9]=1;
bivar[10]=remla[10];bivar[11]=remla[10]*remla[0]/(1-remla[0]);
bivar[12]=remlb[10];bivar[13]=remlb[10]*remlb[0]/(1-remlb[0]);

if(remla[0]>0.001&&remlb[0]>0.001)
{
//gene_bivar(bivar, Ya, Yb, Za, Zb, ZTYa, ZTYb, ZTZa, ZTZb, detZTZa, detZTZb, YTCYa, YTCYb, num_covars, data, respindex, num_samples_use, gmax, keepmults+gstarts[g], weights+gstarts[g]);
}

//print results
fprintf(output3, "%d %s %.3f %.3f %.6f %.6f %.4f %.4e %.4f %.4e %d %.2f\n", g+1, genenames[g], bivar[0], bivar[1], bivar[2], bivar[3], bivar[4], bivar[5], bivar[8], bivar[9], glength, Xsums[num_regs]);
}	//end of mode=11

prev1=gstarts[g];prev2=gends[g];
}	//end of considering gene g
count2++;
}	//end of g loop

fclose(output);if(mode==11){fclose(output2);fclose(output3);}

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}
free(data);free(keepcentres);free(keepmults);

free(Xstarts);free(Xends);free(Xrec);free(Xrev);free(Xsums);free(Xa);
if(num_regs+num_kins==0){free(ZTYa);free(ZTZa);}else{free(Mnull);}
if(mode==11){free(Xb);free(ZTYb);free(ZTZb);}
}

/////////

if(mode==10)	//join gene reml results
{
join_reml(folder, genenames, gstarts, gends, gpartitions, gchr, num_genes, prednames, cut1, cut2, scoretest);
}

if(mode==12)
{
//join bivar
}

/////////

if(mode==16)	//reml analysis
{
if(num_regs>0)
{
Xa=malloc(sizeof(double)*respindex[0][0]*rmax);
Xstarts=malloc(sizeof(int)*num_regs);Xends=malloc(sizeof(int)*num_regs);
Xrec=malloc(sizeof(int)*rmax);Xrev=malloc(sizeof(int)*rmax);
Xsums=malloc(sizeof(double)*num_regs);
fill_X(Xa, Xstarts, Xends, Xrec, Xrev, Xsums, respindex[0], num_samples_use, num_regs, rdata, regindex, rkeepmults, rweights, data, 0, NULL, NULL, &adjust);
}

multi_reml(Ya, Za, Xa, respindex[0][0], num_covars, num_regs, -1, Xstarts, Xends, Xrec, Xrev, Xsums, num_kins, mkins, shortcut, U, E, NULL, outfile, ids1, ids2, rdata_length, rkeepcentres, rkeepmults, rweights, rprednames, ral1, ral2, NULL, NULL, NULL, adjust, adjust2, 0, NULL);
}

if(mode==17)	//decomp kinship and save
{
if((output=fopen(outfile,"w"))==NULL)
{printf("Error writing to %s\n\n", outfile);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output,"%s %s %f ", ids1[i], ids2[i], E[i]);
for(i2=0;i2<num_samples_use;i2++){fprintf(output,"%f ", U[i+i2*num_samples_use]);}
fprintf(output,"\n");
}
fclose(output);
}

if(mode==19||mode==20)	//calc blups or scores
{
//will read in data a bit at a time
data=malloc(sizeof(double)*num_samples_use*bitsize);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

data2=malloc(sizeof(double)*num_samples_use*bitsize);
preds=malloc(sizeof(double*)*(num_kins+num_regs+1));
for(k=0;k<num_kins+num_regs+1;k++)
{
preds[k]=malloc(sizeof(double)*num_samples_use);
for(i=0;i<num_samples_use;i++){preds[k][i]=0;}
}

bittotal=(int)((data_length-1)/bitsize+1);
current=0;weightsum=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(mode==19){printf("Calculating blups for Chunk %d of %d\n", bit+1, bittotal);}
if(mode==20){printf("Calculating scores for Chunk %d of %d\n", bit+1, bittotal);}

current=read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, current, data, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);

for(k=0;k<num_kins+num_regs;k++)	//loop through getting effects and preds
{
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{
if(data[i+j*num_samples_use]!=missingvalue)
{data2[i+j*num_samples_use]=data[i+j*num_samples_use]-allcentres[k][bitstart+j];}
else
{data2[i+j*num_samples_use]=0;}
}
}

if(k<num_kins)	//for kins, get raw effects (must multiply by allfactors twice)
{
alpha=1.0/wsums[k];beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data2, &num_samples_use, mG2[k], &one, &beta, effects[k]+bitstart, &one);
for(j=bitstart;j<bitend;j++){effects[k][j]=effects[k][j]*pow(allfactors[k][j],2);}
}

//now predictions for all
alpha=1.0;beta=1.0;
dgemv_("N", &num_samples_use, &bitlength, &alpha, data2, &num_samples_use, effects[k]+bitstart, &one, &beta, preds[k], &one);
}	//end of k loop
}	//end of bit loop

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}
free(data);free(data2);

for(k=0;k<num_kins+num_regs;k++)	//accumulate (unnecessary for scores)
{
for(i=0;i<num_samples_use;i++){preds[num_kins+num_regs][i]+=preds[k][i];}
for(j=0;j<data_length;j++){effects[num_kins+num_regs][j]+=effects[k][j];}
}

if(mode==19)
{write_blups(effects, data_length, preds, num_samples_use, num_kins, num_regs, prednames, al1, al2, ids1, ids2, allcentres, outfile);}
if(mode==20)
{write_scores(preds, num_samples_use, num_regs, ids1, ids2, outfile);}

for(k=0;k<num_kins+num_regs+1;k++){free(preds[k]);}free(preds);
}	//end of calculating blups

/////////

if(mode==23)
{
manip_kins(outfile, num_kins, kinstems, samedata, kindetails, kinraw, kingz, 1);
}
if(mode==24&&num_regs==0)
{
manip_kins(outfile, num_kins, kinstems, samedata, kindetails, kinraw, kingz, 2);
}
if(mode==24&&num_regs>0)
{
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
read_kinfile(kinstems[0], kins, NULL, num_samples_use, ids1, ids2, 0);

allpreds=malloc(sizeof(char*)*num_predsb);
for(j=0;j<num_predsb;j++){allpreds[j]=malloc(sizeof(char)*100);}
alla1=malloc(sizeof(char)*num_predsb);
alla2=malloc(sizeof(char)*num_predsb);

keepcentres=malloc(sizeof(double)*num_predsb);
keepmults=malloc(sizeof(double)*num_predsb);
weights=malloc(sizeof(float)*num_predsb);
keeppredsb=malloc(sizeof(int)*num_predsb);

weightsum=adjust_details(kinstems[0], allpreds, keepcentres, keepmults, weights, alla1, alla2, keeppredsb, rdata_length, rprednames, rkeepcentres, rkeepmults, rweights, ral1, ral2);

value=weightsum;
for(j=0;j<rdata_length;j++)
{
if(rweights[j]>0)
{
for(i=0;i<num_samples_use;i++)
{
if(rdata[i+j*num_samples_use]!=missingvalue)
{rdata[i+j*num_samples_use]=(rdata[i+j*num_samples_use]-rkeepcentres[j])*rkeepmults[j]*pow(rweights[j],.5);}
else
{rdata[i+j*num_samples_use]=0;}
}
value-=rweights[j];
}
else
for(i=0;i<num_samples_use;i++){rdata[i+j*num_samples_use]=0;}
}

alpha=-1.0/value;beta=weightsum/value;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &rdata_length, &alpha, rdata, &num_samples_use, rdata, &num_samples_use, &beta, kins, &num_samples_use);

write_kins(outfile, kins, NULL, value, ids1, ids2, num_samples_use, keeppredsb, allpreds, keepcentres, keepmults, weights, alla1, alla2, num_predsb, kinraw, kingz);

free(kins);free(alla1);free(alla2);free(keepcentres);free(keepmults);free(weights);free(keeppredsb);
for(j=0;j<num_predsb;j++){free(allpreds[j]);}free(allpreds);
}

/////////

if(mode==25||mode==26||mode==27)	//converting data
{
#include "convert.c"
}	//end of converting data

///////////////////////////

//free variables

//gsl_rng_free(gslr);

if(num_samples!=-1)
{
free(keepsamps);
for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);}
free(ids1);free(ids2);
}

if(num_preds!=-1)
{
free(keeppreds);free(keeppreds_use);
free(chr);free(bp);free(al1);free(al2);
for(j=0;j<data_length;j++){free(prednames[j]);}free(prednames);
}

if(strcmp(weightsfile,"blank")!=0&&(mode==1||mode==2||mode==3))
{free(keeppredsb);}

if(mode==5||mode==7||mode==8||mode==9||mode==11)	
{free(weights);}

if(mode==8||mode==9||mode==10||mode==11||mode==12)
{
for(j=0;j<num_genes;j++){free(genenames[j]);}free(genenames);
free(gchr);
free(gstarts);
free(gends);
free(gpartitions);
}

if(mode==9||mode==11||mode==16||(mode==17&&num_resps>0))
{
free(resp);free(covar);
free(Ya);free(Za);if(mode==11){free(Yb);free(Zb);}
for(m=0;m<num_resps_use;m++){free(respindex[m]);}free(respindex);
}

for(k=0;k<num_kins;k++){free(kinstems[k]);}if(num_kins>0){free(kinstems);}
if(num_kins>0&&(mode==9||mode==16||mode==17))
{for(k=0;k<num_kins;k++){free(mkins[k]);}free(mkins);}

if(mode==17||(shortcut==1&&num_kins==1))
{free(E);free(U);}

if(num_regs>0&&(mode==9||mode==16||mode==24))
{
free(rkeeppreds_use);
for(r=0;r<num_regs;r++){free(regindex[r]);}free(regindex);
for(j=0;j<rdata_length;j++){free(rprednames[j]);}free(rprednames);
free(ral1);free(ral2);
free(rdata);
free(rkeepmults);free(rkeepcentres);free(rweights);
}

if(mode==19||mode==20)
{
for(k=0;k<num_kins+num_regs;k++){free(allcentres[k]);free(allfactors[k]);}
free(allcentres);free(allfactors);free(wsums);
for(k=0;k<num_kins+num_regs+1;k++){free(effects[k]);}free(effects);
for(k=0;k<num_kins;k++){free(mG[k]);free(mG2[k]);}if(num_kins>0){free(mG);free(mG2);}
}

///////////////////////////

printf("Mission completed. All your base are belong to us ;)\n");

return(0);
}	//end of main












