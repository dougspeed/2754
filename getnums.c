/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//see which preds we are using

if(strcmp(mapfile, "blank")!=0){num_preds=countrows(mapfile);}

if(num_preds!=-1)
{
num_preds_use=num_preds;
keeppreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){keeppreds[j]=j;}

if(strcmp(bpredfile,"blank")!=0)
{
allsnps=malloc(sizeof(char*)*num_preds);
for(j=0;j<num_preds;j++){allsnps[j]=malloc(sizeof(char)*100);}
read_strings(mapfile, allsnps, num_preds, 2, 0);

count=countrows(bpredfile);
printf("Reading list of %d predictors to extract from %s\n", count, bpredfile);
usesnps=malloc(sizeof(char*)*count);
for(j=0;j<count;j++){usesnps[j]=malloc(sizeof(char)*100);}
read_strings(bpredfile, usesnps, count, 1, 0);

num_preds_use=find_strings(allsnps, num_preds, usesnps, count, keeppreds, NULL);
if(num_preds_use==0)
{printf("Error, found no predictors in %s in %s\n\n", bpredfile, mapfile);exit(1);}

for(j=0;j<num_preds;j++){free(allsnps[j]);}free(allsnps);
for(j=0;j<count;j++){free(usesnps[j]);}free(usesnps);
}
}	//end of num_preds!=-1


//if re-calculating weights must reduce keeppreds
if(strcmp(weightsfile,"blank")!=0&&(mode==1||mode==2||mode==3))
{
num_preds_useb=num_preds_use;
keeppredsb=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){keeppredsb[j]=keeppreds[j];}

chr=malloc(sizeof(int)*num_preds_useb);
prednames=malloc(sizeof(char*)*num_preds_useb);
for(j=0;j<num_preds_useb;j++){prednames[j]=malloc(sizeof(char)*100);}
bp=malloc(sizeof(double)*num_preds_useb);
al1=malloc(sizeof(char)*num_preds_useb);
al2=malloc(sizeof(char)*num_preds_useb);

if(strcmp(chiamofile,"blank")!=0)
{read_mapfile(mapfile, chr, prednames, bp, al1, al2, num_preds_useb, keeppredsb, 0);}
if(strcmp(chiamofile,"blank")==0)
{read_mapfile(mapfile, chr, prednames, bp, al1, al2, num_preds_useb, keeppredsb, 1);}

weights=malloc(sizeof(double)*num_preds_useb);
read_weightfile(weightsfile, weights, num_preds_useb, keeppredsb, prednames, mapfile);

count=0;
for(j=0;j<num_preds_useb;j++)
{
if(weights[j]>0){keeppreds[count]=keeppredsb[j];count++;}
}
num_preds_use=count;

free(weights);
free(chr);free(bp);free(al1);free(al2);
for(j=0;j<num_preds_useb;j++){free(prednames[j]);}free(prednames);

printf("Re-calcluating weights using the %d (of %d) predictors with non-zero weights from first run\n\n", num_preds_use, num_preds_useb);
}


//and which samples we are using

if(strcmp(famfile,"blank")!=0)
{
num_samples=countrows(famfile);

num_samples_use=num_samples;
keepsamps=malloc(sizeof(int)*num_samples);
for(i=0;i<num_samples;i++){keepsamps[i]=i;}

if(strcmp(bsampfile,"blank")!=0)
{
allsamps=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){allsamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(famfile, allsamps, num_samples);

count=countrows(bsampfile);
printf("Reading list of %d samples to keep from %s\n", count, bsampfile);
usesamps=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){usesamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(bsampfile, usesamps, count);

num_samples_use=find_strings(allsamps, num_samples, usesamps, count, keepsamps, NULL);
if(num_samples_use==0)
{printf("Error, found no samples in %s in %s\n\n", bsampfile, famfile);exit(1);}

for(i=0;i<num_samples;i++){free(allsamps[i]);}free(allsamps);
for(i=0;i<count;i++){free(usesamps[i]);}free(usesamps);
}

//if subtracting regions, make sure individuals match first kinship file ids
if(mode==24&&num_regs>0)
{
allsamps=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){allsamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(famfile, allsamps, num_samples);

sprintf(filename,"%s.grm.id", kinstems[0]);
count=countrows(filename);
usesamps=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){usesamps[i]=malloc(sizeof(char)*200);}
read_ids_merge(filename, usesamps, count);

num_samples_use=find_strings(allsamps, num_samples, usesamps, count, keepsamps, NULL);
if(num_samples_use!=count)
{printf("Error, only %d of the %d samples in %s are in %s\n\n", num_samples_use, count, filename, famfile);exit(1);}

for(i=0;i<num_samples;i++){free(allsamps[i]);}free(allsamps);
for(i=0;i<count;i++){free(usesamps[i]);}free(usesamps);
}	
}




