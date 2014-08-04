/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//user-defined variables

int mode=0;

char folder3[500]="blank", folder2[500]="", folder[500]="";
char outfile2[500]="", outfile[500]="blank";

char workdir2[500]="blank", workdir[500]="";

char bedfile2[500]="blank", bedfile[500]="";
char chiamofile3[500]="blank", chiamofile2[500]="", chiamofile[500]="";
char spfile2[500]="blank", spfile[500]="";
char speedfile2[500]="blank", speedfile[500]="";
char mapfile2[500]="blank", mapfile[500]="";
char famfile2[500]="blank", famfile[500]="";
char datafile3[500]="blank", datafile2[500]="blank", datafile[500]="";

int chiamoprobs=3;
int chiamoheaders=5;
char chiamosuffix[500]="gen";
double missingvalue=-99;
float power=-1;

//char kpredfile2[500]="blank", kpredfile[500]="";
char bpredfile2[500]="blank", bpredfile[500]="";
//char ksampfile2[500]="blank", ksampfile[500]="";
char bsampfile2[500]="blank", bsampfile[500]="";

int num_subs=-1;
char subpref2[500]="blank", subpref[500];

int num_parts=0;
char partpref2[500]="blank", partpref[500];

int num_regs=0;
char regpref2[500]="blank", regpref[500];

float minmaf=-1;
float maxmaf=-1;
float minvar=-1;
float minobs=-1;
int encoding=1;
float threshold=-1;

char weightsfile2[500]="blank", weightsfile[500]="";
int ignoreweights=0;

int section=-1;
int section_length=3000;
int section_buffer=500;

int partition=-1;
int partition_length=1000000;
int bychr=0;
int kinraw=0;
int kingz=0;

char genefile2[500]="blank", genefile[500]="";
int chunksbp=-1;
double chunks=-1;
int gene_buffer=0;
float minweight=0;
int overlap=0;
char pvafile2[500]="blank", pvafile[500]="";

char respfile2[500]="blank", respfile[500]="";
int mpheno=-1, mpheno2=-1;
int permute=0;

char covarfile2[500]="blank", covarfile[500]="";

char kinname2[500]="blank", kinname[500]="";
char kinlist2[500]="blank", kinlist[500]="";
char eigenfile2[500]="blank", eigenfile[500]="";

char remlfile2[500]="blank", remlfile[500]="";
char blupfile2[500]="blank", blupfile[500]="";
char regfile2[500]="blank", regfile[500]="";

int num_snps=-1;
int num_phenos=100;
int num_causals=100;
float her=-1;
char causalsfile[500], causalsfile2[500]="blank";
char effectsfile[500], effectsfile2[500]="blank";

int num_kins=0;
char **kinstems;
int kindetails=1;

int ignoremissing=1;	//I might fix this to 1

float priora=-1.0, priorb=-1.0;	//default is -1 and 1
float adjust=0.0;

float cut1=-1, cut2=-1;
int scoretest=0;


///////////////////////////

//program variables

int i, i2, j, j2, k, g, m, r;
int count, found, found2, found3, current, prev1, prev2, length, lwork, info, one=1;
int readint, *order, gen, value, conv[4], shortcut=0;
double alpha, beta, wkopt, *work;

char readstring[100], readstring2[100], readstring3[100];

char filename[500], filename2[500], filename3[500];
char argname[50];

FILE *input, *datainput, *output, *output2, *output3;

int evalue, evalue2;	//these are used for checking folder exists
struct stat statstruct;


int num_samples=-1, num_samples_use=-1;
int num_preds=-1, num_preds_use=-1, num_preds_useb;

char **allsnps, **usesnps;
char **allsamps, **usesamps;

int *keeppreds, *keeppreds_use, *keeppredsb, *rkeeppreds_use;
int *keepsamps;

char **ids1, **ids2;

int *chr, *rchr;
char **prednames, **rprednames, *al1, *al2, *ral1, *ral2;
double *bp, *rbp;

float *weights, *rweights, weightsum;

int data_length, rdata_length;
double *data, *rdata, *data2;
double *keepcentres, *keepmults, *rkeepcentres, *rkeepmults;

double *kins, *divs, **mkins, *U, *E;

int bit, bittotal, bitstart, bitend, bitlength, bitsize=5000;

int num_genes, num_genes_approx, gene_max=0;
int *gchr, *gstarts, *gends, *gpartitions;
double *gbp1, *gbp2;
char **genenames;

int num_resps, num_resps_use, keepresps[2], **respindex;
double *resp, *resptemp;

int **regindex;

int num_covars=1;
double *covar;

double *Ya, *Za, *ZTYa, *ZTZa, YTYa, detZTZa, YTCYa, remla[12];
double *Yb, *Zb, *ZTYb, *ZTZb, YTYb, detZTZb, YTCYb, remlb[12];
double bivar[11];

double **mG, **mG2, **effects, **preds;
float **allcentres, **allfactors, *wsums;




