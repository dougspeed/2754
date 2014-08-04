/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int set_up_sections(char *folder, int *chr, double *bp, int num_preds_use, int section_length, int buffer, char *mapfile, char *datafile, char *kpredfile, int type)
{
int j;
int numchr, *chrlengths, *chrstarts;
int start1, end1, start2, end2;
int section, num_sections, length_sections, total_sections;

char filename[500];
FILE *output;


//find total number of chr
numchr=0;
for(j=0;j<num_preds_use;j++){if(chr[j]>numchr){numchr=chr[j];}}
numchr++;

//find numbers of predictors on each chr and starts
chrlengths=malloc(sizeof(int)*numchr);
chrstarts=malloc(sizeof(int)*numchr);
for(j=0;j<numchr;j++){chrlengths[j]=0;}
for(j=num_preds_use-1;j>=0;j--){chrlengths[chr[j]]++;chrstarts[chr[j]]=j;}

//open file which stores section details and write header lines
sprintf(filename,"%ssection_details.txt",folder);
if(type==1){sprintf(filename,"%sre-section_details.txt",folder);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s - check write permission granted\n\n",filename);exit(1);}

fprintf(output, "Datafiles: %s %s\n", datafile, mapfile);

if(strcmp(kpredfile,"blank")==0){fprintf(output, "Using All Predictors\n");}
else{fprintf(output, "Using keepfile %s\n", kpredfile);}
fprintf(output,"Section Analysed_Predictors Retained_Predictors\n");

//loop through chromosomes
total_sections=0;
for(j=0;j<numchr;j++)
{
if(chrlengths[j]>0)
{
//get the lengths
num_sections=(chrlengths[j]-1)/section_length+1;
if(num_preds_use<=section_length+2*buffer){num_sections=1;}
length_sections=(chrlengths[j]-1)/num_sections+1;

for(section=0;section<num_sections;section++)
{
//section will contain the boundaries of 1:length_sections predictors + buffer predictors, and the central ones
start2=length_sections*section;
end2=length_sections*(section+1);
start1=start2-buffer;
end1=end2+buffer;

//check all these are within 0:(num_preds_use-1);
if(start1<0){start1=0;}
if(end1>chrlengths[j]){end1=chrlengths[j];}
if(end2>chrlengths[j]){end2=chrlengths[j];}

//write predictors to analyse, then to retain
fprintf(output,"%d %d-%d %d-%d\n", total_sections+1, chrstarts[j]+start1+1, chrstarts[j]+end1, start2-start1+1, end2-start1);

total_sections++;
}	//end of sections loop
}
}

fclose(output);

printf("The %d predictors to be examined have been split into %d sections of (approx) length %d\n", num_preds_use, total_sections, section_length);
printf("The details for sections are stored in %s\n\n", filename);

if(strcmp(kpredfile,"blank")!=0)
{printf("Weightings will be calculated only for subset S of all predictors\nPredictors not in S can not be used when calculating kinships\nIf you later wish to calculate kinships over subsets, make sure these subsets are contained within S.\n\n");}

free(chrlengths);free(chrstarts);

return(0);
}	//end of set_up_sections


///////////////////////////


int get_section_boundaries(char *folder, int section, int *keeppreds_use, int *keeppreds, int num_preds_use, int type)
{
int j;
int section2, start, end;
int ignore;

char readchar;

char filename[500];
FILE *input;


//first set the details filename
sprintf(filename,"%ssection_details.txt", folder);
if(type==1){sprintf(filename,"%sre-section_details.txt", folder);}

//check section is reasonable
if(section>countrows(filename)-3)
{printf("The section number specified (%d) is larger than the total number of sections (%d) listed in %s\n\n", section, countrows(filename)-3, filename);exit(1);}

//open the file
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

//to find the start and end for section, skip first sections+2 rows
section2=0;
while(section2<section+2)
{
readchar=0;
while(readchar!=10){(void)fscanf(input,"%c", &readchar);}
section2++;
}

//read the row
if(fscanf(input,"%d %d-%d %d-%d ", &section2, &start, &end, &ignore, &ignore)!=5)
{printf("Error reading details for Section %d (Row %d) from %s\n\n", section, section+3, filename);exit(1);}
if(section2!=section){printf("Error reading %s; Row %d should contain details for section %d (not %d)\n\n", filename, section+3, section, section2);exit(1);}
start=start-1;

//check sensible
if(end>num_preds_use)
{printf("Error reading %s; the end predictor (%d) of Section %d is greater than the total number of predictors being considered (%d). Has file been altered since creation using argument \"--cut-weights\"\n\n", filename, end, section, num_preds_use);exit(1);}

fclose(input);

//load up keeppreds_use
for(j=0;j<end-start;j++){keeppreds_use[j]=keeppreds[start+j];}

return(end-start);
}	//end of get_section_boundaries


///////////////////////////



///////////////////////////





