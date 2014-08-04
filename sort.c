/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//these are used with qsort???

int compare_int (const void *a, const void *b)
{
return(*(int*)a-*(int*)b);
}

int compare_float (const void *a, const void *b)
{
return(*(float*)a-*(float*)b);
}

int compare_double (const void *a, const void *b)
{
return(*(double*)a-*(double*)b);
}

int compare_string (const void *a, const void *b)
{
return(strcmp(*(char**)a,*(char**)b));
}

/////////

void swap_int(int *a, int *b)
{
int c;c=*a;*a=*b;*b=c;
}

void swap_double(double *a, double *b)
{
double c;c=*a;*a=*b;*b=c;
}

/////////

void pivotindex_int(int *arr, int *indindex, int beg, int end)
{
int l, r;
l=beg+1; r=end+1;
while(l<r)
{
if(arr[indindex[l]]<arr[indindex[beg]]){l++;}
else{r--;swap_int(indindex+l, indindex+r);}
}
r--;
swap_int(indindex+beg, indindex+r);
if(beg<r-1){pivotindex_int(arr, indindex, beg, r-1);}
if(r+1<end){pivotindex_int(arr, indindex, r+1, end);}
}

void sortindex_int(int *arr, int *indindex, int length)
{
int i;
for(i=0;i<length;i++){indindex[i]=i;}
pivotindex_int(arr, indindex, 0, length-1);
}

/////////

void pivotindex_string(char **arr, int *indindex, int beg, int end)
{
int l, r;
l=beg+1; r=end+1;
while(l<r)
{
if(strcmp(arr[indindex[l]],arr[indindex[beg]])<=0){l++;}
else{r--;swap_int(indindex+l, indindex+r);}
}
r--;
swap_int(indindex+beg, indindex+r);
if(beg<r-1){pivotindex_string(arr, indindex, beg, r-1);}
if(r+1<end){pivotindex_string(arr, indindex, r+1, end);}
}

void sortindex_string(char **arr, int *indindex, int length)
{
int i;
for(i=0;i<length;i++){indindex[i]=i;}
pivotindex_string(arr, indindex, 0, length-1);
}

/////////

int find_strings(char **str1, int length1, char **str2, int length2, int *keepindex, int *keepindexb)
//sees which of str1 are in str2 - gets index of str1 used and to which these match in str2
{
int j, j2, count;
int *index1, *index2, *found;

//first sort the two strings
index1=malloc(sizeof(int)*length1);
index2=malloc(sizeof(int)*length2);

sortindex_string(str1, index1, length1);
sortindex_string(str2, index2, length2);

found=malloc(sizeof(int)*length1);
count=0;j2=0;
for(j=0;j<length1;j++)
{
found[index1[j]]=-1;
while(j2<length2)
{
if(strcmp(str1[index1[j]],str2[index2[j2]])==0)
{found[index1[j]]=index2[j2];j2++;break;}
if(strcmp(str1[index1[j]],str2[index2[j2]])<0)
{break;}
j2++;
}
}

count=0;
for(j=0;j<length1;j++)
{
if(found[j]!=-1)
{
keepindex[count]=j;
if(keepindexb!=NULL){keepindexb[count]=found[j];}
count++;
}
}

free(index1);free(index2);free(found);

return(count);
}

/////////

int find_ids(int *keepindex, char **ids1, char **ids2, int num_samples_use, char **ids1b, char **ids2b, int length)
{
int i, ib, found;

for(i=0;i<num_samples_use;i++)
{
found=0;
for(ib=0;ib<length;ib++)
{
if(strcmp(ids1[i],ids1b[ib])==0&&strcmp(ids2[i],ids2b[ib])==0)
{
if(found!=0)	//already found
{return(-(i+1));}
else
{keepindex[i]=ib;found=1;}
}
}
if(found==0)	//can't find
{return(i+1);}
}

return(0);
}

/////////

int update_resps(double *resp, int **respindex, char **ids1, char **ids2, int *keepsamps, int num_samples_use, int num_resps_use)
{
int i, m, count;
int *found, *rev, *kindex;
double *respb;
char **ids1b, **ids2b;

//any need to do this?
if(respindex[0][0]==num_samples_use){return(num_samples_use);}
if(num_resps_use==2)
{
if(respindex[1][0]==num_samples_use){return(num_samples_use);}
}

//find union of all samples
found=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){found[i]=0;}
kindex=malloc(sizeof(int)*num_samples_use);
rev=malloc(sizeof(int)*num_samples_use);

for(m=0;m<num_resps_use;m++)
{
for(i=0;i<respindex[m][0];i++){found[respindex[m][1+i]]++;}
}
count=0;
for(i=0;i<num_samples_use;i++)
{
if(found[i]>0){kindex[count]=i;rev[i]=count;count++;}
}

if(count<num_samples_use)
{
respb=malloc(sizeof(double)*num_samples_use*num_resps_use);
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){respb[i+m*num_samples_use]=resp[i+m*num_samples_use];}
}
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<count;i++){resp[i+m*count]=respb[kindex[i]+m*num_samples_use];}
}
free(respb);

for(m=0;m<num_resps_use;m++)
{
for(i=0;i<respindex[m][0];i++){respindex[m][1+i]=rev[respindex[m][1+i]];}
}

ids1b=malloc(sizeof(char*)*num_samples_use);
ids2b=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{ids1b[i]=malloc(sizeof(char)*100);strcpy(ids1b[i],ids1[i]);
ids2b[i]=malloc(sizeof(char)*100);strcpy(ids2b[i],ids2[i]);}

for(i=0;i<count;i++){strcpy(ids1[i],ids1b[kindex[i]]);strcpy(ids2[i],ids2b[kindex[i]]);}
for(i=count;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);}
for(i=0;i<num_samples_use;i++){free(ids1b[i]);free(ids2b[i]);}
free(ids1b);free(ids2b);

for(i=0;i<count;i++){keepsamps[i]=keepsamps[kindex[i]];}

printf("Due to missing phenotypic values, number of samples reduced to %d\n\n", count);
}	//end of reducing

free(found);free(rev);free(kindex);

return(count);
}	//end of update_resps








