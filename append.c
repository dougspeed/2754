/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//complete suffix of chiamofile (and datafile) and provide a warning
if(strcmp(chiamofile3,"blank")!=0)
{
sprintf(chiamofile2,"%s.%s", chiamofile3, chiamosuffix);
sprintf(datafile2,"%s.%s", datafile3, chiamosuffix);
printf("When using a chiamo file, the default row format is 5 header columns\n(2xSNP_IDs, BP, 2xAlleles) and 3 probabilities per sample\nThese can be changed with the arguments \"--chiamo-headers\" and \"--chiamo-probs\"\nThe default extension is \".gen\" which can be changed using \"--chiamo-suffix\"\n\n");
}
else{strcpy(chiamofile2,"blank");}


//if workdir has been specified, make sure it ends in "/"
if(strcmp(workdir2,"blank")!=0)
{
count=strlen(workdir2);
if(workdir2[count-1]!='/')
{sprintf(workdir,"%s/",workdir2);}
else
{strcpy(workdir,workdir2);}
}
else
{
strcpy(workdir,"");
}


//make sure folder ends in "/"
if(strcmp(folder3,"blank")!=0)
{
count=strlen(folder3);
if(folder3[count-1]!='/')
{sprintf(folder2,"%s/",folder3);}
else
{strcpy(folder2,folder3);}
}


//prefix folder and files with working directory, unless they are absolute filenames

if(strcmp(folder3,"blank")!=0)
{
if(folder2[0]!='/')
{sprintf(folder,"%s%s", workdir, folder2);}
else
{strcpy(folder,folder2);}
}

if(outfile2[0]!='/')
{sprintf(outfile,"%s%s", workdir, outfile2);}
else
{strcpy(outfile,outfile2);}

/////////

if(strcmp(bedfile2,"blank")!=0)
{
if(bedfile2[0]!='/')
{sprintf(bedfile,"%s%s", workdir, bedfile2);}
else
{strcpy(bedfile,bedfile2);}
}
else{strcpy(bedfile,"blank");}

if(strcmp(chiamofile2,"blank")!=0)
{
if(chiamofile2[0]!='/')
{sprintf(chiamofile,"%s%s", workdir, chiamofile2);}
else
{strcpy(chiamofile,chiamofile2);}
}
else{strcpy(chiamofile,"blank");}

if(strcmp(spfile2,"blank")!=0)
{
if(spfile2[0]!='/')
{sprintf(spfile,"%s%s", workdir, spfile2);}
else
{strcpy(spfile,spfile2);}
}
else{strcpy(spfile,"blank");}

if(strcmp(speedfile2,"blank")!=0)
{
if(speedfile2[0]!='/')
{sprintf(speedfile,"%s%s", workdir, speedfile2);}
else
{strcpy(speedfile,speedfile2);}
}
else{strcpy(speedfile,"blank");}

/////////

if(strcmp(mapfile2,"blank")!=0)
{
if(mapfile2[0]!='/')
{sprintf(mapfile,"%s%s", workdir, mapfile2);}
else
{strcpy(mapfile,mapfile2);}
}
else{strcpy(mapfile,"blank");}

if(strcmp(famfile2,"blank")!=0)
{
if(famfile2[0]!='/')
{sprintf(famfile,"%s%s", workdir, famfile2);}
else
{strcpy(famfile,famfile2);}
}
else{strcpy(famfile,"blank");}

if(strcmp(datafile2,"blank")!=0)
{
if(datafile2[0]!='/')
{sprintf(datafile,"%s%s", workdir, datafile2);}
else
{strcpy(datafile,datafile2);}
}
else{strcpy(datafile,"blank");}

/////////

if(strcmp(weightsfile2,"blank")!=0)
{
if(weightsfile2[0]!='/')
{sprintf(weightsfile,"%s%s", workdir, weightsfile2);}
else
{strcpy(weightsfile,weightsfile2);}
}
else{strcpy(weightsfile,"blank");}

if(strcmp(respfile2,"blank")!=0)
{
if(respfile2[0]!='/')
{sprintf(respfile,"%s%s", workdir, respfile2);}
else
{strcpy(respfile,respfile2);}
}
else{strcpy(respfile,"blank");}

if(strcmp(covarfile2,"blank")!=0)
{
if(covarfile2[0]!='/')
{sprintf(covarfile,"%s%s", workdir, covarfile2);}
else
{strcpy(covarfile,covarfile2);}
}
else{strcpy(covarfile,"blank");}

if(strcmp(genefile2,"blank")!=0)
{
if(genefile2[0]!='/')
{sprintf(genefile,"%s%s", workdir, genefile2);}
else
{strcpy(genefile,genefile2);}
}
else{strcpy(genefile,"blank");}

if(strcmp(pvafile2,"blank")!=0)
{
if(pvafile2[0]!='/')
{sprintf(pvafile,"%s%s", workdir, pvafile2);}
else
{strcpy(pvafile,pvafile2);}
}
else{strcpy(pvafile,"blank");}

/////////
/*
if(strcmp(kpredfile2,"blank")!=0)
{
if(kpredfile2[0]!='/')
{sprintf(kpredfile,"%s%s", workdir, kpredfile2);}
else
{strcpy(kpredfile,kpredfile2);}
}
else{strcpy(kpredfile,"blank");}

if(strcmp(ksampfile2,"blank")!=0)
{
if(ksampfile2[0]!='/')
{sprintf(ksampfile,"%s%s", workdir, ksampfile2);}
else
{strcpy(ksampfile,ksampfile2);}
}
else{strcpy(ksampfile,"blank");}
*/

if(strcmp(bpredfile2,"blank")!=0)
{
if(bpredfile2[0]!='/')
{sprintf(bpredfile,"%s%s", workdir, bpredfile2);}
else
{strcpy(bpredfile,bpredfile2);}
}
else{strcpy(bpredfile,"blank");}

if(strcmp(bsampfile2,"blank")!=0)
{
if(bsampfile2[0]!='/')
{sprintf(bsampfile,"%s%s", workdir, bsampfile2);}
else
{strcpy(bsampfile,bsampfile2);}
}
else{strcpy(bsampfile,"blank");}

/////////

if(strcmp(kinname2,"blank")!=0)
{
if(kinname2[0]!='/')
{sprintf(kinname,"%s%s", workdir, kinname2);}
else
{strcpy(kinname,kinname2);}
}
else{strcpy(kinname,"blank");}

if(strcmp(kinlist2,"blank")!=0)
{
if(kinlist2[0]!='/')
{sprintf(kinlist,"%s%s", workdir, kinlist2);}
else
{strcpy(kinlist,kinlist2);}
}
else{strcpy(kinlist,"blank");}

if(strcmp(eigenfile2,"blank")!=0)
{
if(eigenfile2[0]!='/')
{sprintf(eigenfile,"%s%s", workdir, eigenfile2);}
else
{strcpy(eigenfile,eigenfile2);}
}
else{strcpy(eigenfile,"blank");}

if(strcmp(remlfile2,"blank")!=0)
{
if(remlfile2[0]!='/')
{sprintf(remlfile,"%s%s", workdir, remlfile2);}
else
{strcpy(remlfile,remlfile2);}
}
else{strcpy(remlfile,"blank");}

if(strcmp(blupfile2,"blank")!=0)
{
if(blupfile2[0]!='/')
{sprintf(blupfile,"%s%s", workdir, blupfile2);}
else
{strcpy(blupfile,blupfile2);}
}
else{strcpy(blupfile,"blank");}

if(strcmp(regfile2,"blank")!=0)
{
if(regfile2[0]!='/')
{sprintf(regfile,"%s%s", workdir, regfile2);}
else
{strcpy(regfile,regfile2);}
}
else{strcpy(regfile,"blank");}

/////////

if(strcmp(subpref2,"blank")!=0)
{
if(subpref2[0]!='/')
{sprintf(subpref,"%s%s", workdir, subpref2);}
else
{strcpy(subpref,subpref2);}
}
else{strcpy(subpref,"blank");}

if(strcmp(partpref2,"blank")!=0)
{
if(partpref2[0]!='/')
{sprintf(partpref,"%s%s", workdir, partpref2);}
else
{strcpy(partpref,partpref2);}
}
else{strcpy(partpref,"blank");}

if(strcmp(regpref2,"blank")!=0)
{
if(regpref2[0]!='/')
{sprintf(regpref,"%s%s", workdir, regpref2);}
else
{strcpy(regpref,regpref2);}
}
else{strcpy(regpref,"blank");}




