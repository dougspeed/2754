

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<errno.h>
#include<time.h>


//use one of these two includes if you have glpk installed on your computer
#include<glpk.h>
//#include<glpk/glpk.h>

//and use these includes if you have gsl installed on your computer
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


//will export functions from c

//extern void dgemm_();
//extern void dgemv_();
//extern void dpotrf_();
//extern void dpotri_();
//extern void dpotrs_();
//extern void dgetrf_();
//extern void dgetrs_();
//extern void dsyev_();


#include "sort.c"
#include "norm.c"
#include "fileops.c"



int main (int argc, const char * argv[])
{
int i, j;
int start, end;

int num_samples, num_samples_use, num_preds, data_length;
int *keepsamps, *keeppreds;

double *dataraw, *datastan;

char bedfile[500], bimfile[500], famfile[500];

strcpy(bedfile,"data.bed");
strcpy(bimfile,"data.bim");
strcpy(famfile,"data.fam");

start=0;end=100;

num_samples=countrows(famfile);
num_preds=countrows(bimfile);

num_samples_use=500;data_length=100;

keepsamps=malloc(sizeof(int)*num_samples_use);
for(i=0;i<500;i++){keepsamps[i]=4*i;}

keeppreds=malloc(sizeof(int)*data_length);
for(j=0;j<data_length;j++){keeppreds[j]=j*2+1;}

dataraw=malloc(sizeof(double)*num_samples_use*data_length);

if(strcmp("bedfile","blank")!=0)
{read_bed_fly(bedfile, dataraw, num_samples_use, keepsamps, start, end, keeppreds, num_samples, num_preds, -99.0);}



}









