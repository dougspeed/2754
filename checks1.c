/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//check that specified folder exists
if(strcmp(folder3,"blank")!=0)
{
evalue=stat(folder, &statstruct);
evalue2=0;
if(evalue!=0)
{evalue2=mkdir(folder, S_IRWXU | S_IROTH);}
if(evalue2!=0){printf("Error, failed to create folder %s, check write permissions or that a file of that name does not already exist\n\n",folder);exit(1);}
}


if(mode==2||mode==3)
{printf("Not possible to calculate weights using this version of LDAK\n\n");exit(1);}


//checks for data - don't need for modes 3, 6, 10+nocut, 12, 16+noregs, 17, 18, 22, 23 and 24+noregs
//although can do without datafiles for modes 1, 4 and 7

if(found2>1)	//too many data arguments provided
{printf("Error, can not provide more than one argument from \"--bfile\", \"--chiamo\", \"--sp\" or \"--speed\"\n\n");exit(1);}

//check datafiles provided when required
if(mode==1||mode==2||mode==4||mode==5||mode==7||mode==8||mode==9||(mode==10&&cut1!=-1)||mode==11||(mode==16&&num_regs>0)||mode==19||mode==20||mode==21||(mode==24&&num_regs>0)||mode==25||mode==26||mode==27)
{
if(found2==0)
{printf("Error, necessary to provide prefix for datafiles using one argument from \"--bfile\", \"--chiamo\", \"--sp\" or \"--speed\"\n\n");exit(1);}
}

//check mapfile provided and exists when required
if(mode==1||mode==2||mode==4||mode==5||mode==7||mode==8||mode==9||(mode==10&&cut1!=-1)||mode==11||(mode==16&&num_regs>0)||mode==19||mode==20||mode==21||(mode==24&&num_regs>0)||mode==25||mode==26||mode==27)
{
if((input=fopen(mapfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",mapfile);exit(1);}
fclose(input);
}

//check bedfile/chiamofile/spfile/speedfile and famfile exist
if(mode==2||mode==5||mode==8||mode==9||mode==11||(mode==16&&num_regs>0)||mode==19||mode==20||mode==21||(mode==24&&num_regs>0)||mode==25||mode==26||mode==27)
{
if(strcmp(bedfile,"blank")!=0)
{
if((input=fopen(bedfile,"rb"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",bedfile);exit(1);}
fclose(input);
}
if(strcmp(chiamofile,"blank")!=0)
{
if((input=fopen(chiamofile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",chiamofile);exit(1);}
fclose(input);
}
if(strcmp(spfile,"blank")!=0)
{
if((input=fopen(spfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",spfile);exit(1);}
fclose(input);
}
if(strcmp(speedfile,"blank")!=0)
{
if((input=fopen(speedfile,"rb"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",speedfile);exit(1);}
fclose(input);
}

if((input=fopen(famfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",famfile);exit(1);}
fclose(input);
}

///////////////////////////

//check for keepfiles

if(strcmp(bsampfile,"blank")!=0)	//check bsampfile exists
{
if((input=fopen(bsampfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",bsampfile);exit(1);}
fclose(input);
}

if(strcmp(bpredfile,"blank")!=0)	//check bpredfile exists
{
if((input=fopen(bpredfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",bpredfile);exit(1);}
fclose(input);
}

///////////////////////////

//check for kinships and sort

if(strcmp(kinname,"blank")!=0&&strcmp(kinlist,"blank")!=0)
{printf("Error, can not supply both a kinfile (\"--grm\") and a kinlist (\"--mgrm\")\n\n");exit(1);}

//won't check whether kinships can (will) be used

if(strcmp(kinname,"blank")!=0)	//will have already been prefixed by folder
{
num_kins=1;
kinstems=malloc(sizeof(char *)*num_kins);
kinstems[0]=malloc(sizeof(char)*500);
strcpy(kinstems[0], kinname);
}

if(strcmp(kinlist,"blank")!=0)
{
if((input=fopen(kinlist,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",kinlist);exit(1);}
fclose(input);

num_kins=countrows(kinlist);
kinstems=malloc(sizeof(char *)*num_kins);

if((input=fopen(kinlist,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",kinlist);exit(1);}
for(k=0;k<num_kins;k++)
{
kinstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, kinlist);exit(1);}

if(readstring[0]!='/'){sprintf(kinstems[k], "%s%s", workdir, readstring);}
else{strcpy(kinstems[k], readstring);}
}
fclose(input);
}

if(mode==6)	//joining kins
{
sprintf(filename, "%spartition_details.txt", folder);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",filename);exit(1);}

num_kins=countrows(filename)-3;
kinstems=malloc(sizeof(char *)*num_kins);

for(k=0;k<num_kins;k++)
{kinstems[k]=malloc(sizeof(char)*500);sprintf(kinstems[k], "%skinship%d", folder, k+1);}
sprintf(outfile, "%skinshipALL", folder);
}

//now check present
for(k=0;k<num_kins;k++)
{
sprintf(filename,"%s.grm.id", kinstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",filename);exit(1);}
fclose(input);

sprintf(filename,"%s.grm.bin", kinstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",filename);exit(1);}
fclose(input);

sprintf(filename,"%s.grm.N.bin", kinstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",filename);exit(1);}
fclose(input);

if(kindetails==1)
{
sprintf(filename,"%s.grm.details", kinstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, %s does not exist - should have been created using argument \"--calc-kins\"\nTo proceed regardless, use option \"--kinship-details NO\"\n\n",filename);exit(1);}
fclose(input);

count=countcols(filename);
if(count!=7)
{printf("Error, %s should have 7 columns (not %d)\n\n", filename, count);exit(1);}

sprintf(filename,"%s.grm.adjust", kinstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, %s does not exist - should have been created using argument \"--calc-kins\"\nTo proceed regardless, use option \"--no-details YES\"\n\n",filename);exit(1);}
fclose(input);
}
}

if(strcmp(eigenfile,"blank")!=0)
{
if(mode!=16)
{printf("Error, can only provide eigen-decomposition with argument \"--reml\"\n\n");exit(1);}

if(num_kins!=1)
{printf("Error, can only provide eigen-decomposition when using one kinship (not %d)\n\n", num_kins);exit(1);}

if((input=fopen(eigenfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",eigenfile);exit(1);}
fclose(input);
}

///////////////////////////

//check for subset files

if(num_subs!=-1&&strcmp(subpref,"blank")==0)
{printf("If using the argument \"--subset-number\", you must also use \"--subset-prefix\"\n\n");exit(1);}

if(num_subs==-1&&strcmp(subpref,"blank")!=0)
{printf("If using the argument \"--subset-prefix\", you must also use \"--subset-number\"\n\n");exit(1);}

if(num_subs!=-1)		//check allowed subsets, and subpref provided
{
if(mode!=1&&mode!=2&&mode!=3)
{printf("Error, only allowed to provide subsets when dealing with weights\n\n");exit(1);}

printf("Will get subsets of individuals (e.g. cases and controls) from %s%d, ... %s%d\n", subpref, 1, subpref, num_subs);
printf("Note that now, these files must contain IDs (not indexes). Sorry!\n\n");exit(1);
}
else
{
num_subs=1;
if(mode==1||mode==2||mode==3)
{printf("If using a binary outcome and samples genotyped separately, advisable to use Subset Options\nSee http://dougspeed.com/ldak/subsets for details.\n\n");}
}

///////////////////////////

//check for partition files and sort

if(num_parts!=0&&strcmp(partpref,"blank")==0)
{printf("If using the argument \"--partition-number\", you must also use \"--partition-prefix\"\n\n");exit(1);}

if(num_parts==0&&strcmp(partpref,"blank")!=0)
{printf("If using the argument \"--partition-prefix\", you must also use \"--partition-number\"\n\n");exit(1);}

if(num_parts>0)
{
if(mode!=4&&mode!=5&&mode!=6)
{printf("Error, only allowed to provide partitions when calculating kinships\n\n");exit(1);}

if(strcmp(bpredfile,"blank")!=0)
{printf("Warning, can not extract predictors if providing partitions\n\n");exit(1);}

for(j=0;j<num_parts;j++)
{
sprintf(filename, "%s%d", partpref, j+1);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, file containing predictors for Partition %d (%s) does not exist.\n\n", j+1, filename);exit(1);}
fclose(input);

if(countrows(filename)<1)
{printf("Error, file containing predictors for Partition %d (%s) contains no predictors.\n\n", j+1, filename);exit(1);}
}
}

///////////////////////////

//check for region files

if(num_regs!=0&&strcmp(regpref,"blank")==0)
{printf("If using the argument \"--region-number\", you must also use \"--region-prefix\"\n\n");exit(1);}

//if(num_regs==0&&strcmp(regpref,"blank")!=0)
//{printf("If using the argument \"--region-prefix\", you must also use \"--region-number\"\n\n");exit(1);}

if(num_regs>0)
{
if(mode==19)
{printf("Error, when calculating blups with regions use \"--regfile\" instead of \"--region-number\" and \"--region-prefix\"\n\n");exit(1);}

if(mode!=9&&mode!=16&&mode!=24)
{printf("Error, can only provide region predictors with \"--calc-genes-reml\", \"--reml\" or \"--sub-grm\"\n\n");exit(1);}

if(strcmp(bpredfile,"blank")!=0)
{printf("Error, can not extract predictors if providing regions\n\n");exit(1);}

for(r=0;r<num_regs;r++)
{
sprintf(filename, "%s%d", regpref, r+1);
if((input=fopen(filename,"r"))==NULL)
{printf("Error, file containing predictors for Region %d (%s) does not exist.\n\n", r+1, filename);exit(1);}
fclose(input);

if(countrows(filename)<1)
{printf("Error, file containing predictors for Region %d (%s) contains no predictors.\n\n", r+1, filename);exit(1);}
}
}

///////////////////////////

//checks for weightfile

if(strcmp(weightsfile,"blank")!=0&&ignoreweights==1)
{printf("It is contradictory to both provide a weights file using argument \"--weights\" and use the argument \"--ignore-weights YES\"\n\n");exit(1);}

//won't check whether weights can (will) be used

if(strcmp(weightsfile,"blank")!=0)	//check weightsfile exists and right size
{
if((input=fopen(weightsfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",weightsfile);exit(1);}
fclose(input);

if(found2==0)
{printf("Error; to use weights must also provide prefix for datafiles\n\n");exit(1);}

if(countrows(weightsfile)!=countrows(mapfile))
{printf("Error; number of rows of %s (%d) should equal number of rows of %s (%d)\n\n", weightsfile, countrows(weightsfile), mapfile, countrows(mapfile));exit(1);}

count=countcols(weightsfile);
if(count==4)
{printf("Error, it looks like you have provided a weights file from an old version of LDAK; now the file must have either 5 columns (with predictor names in final column) or just 1 column. Sorry!\n\n");exit(1);}
if(count!=1&&count!=5)
{printf("Error, %s has %d columns\nThe default format (from \"--join-weights\") is 5 columns; although a single column file can be used\n\n", weightsfile, count);exit(1);}
if(count==1)
{printf("Warning, %s contains only weights (not predictor names) so make sure these align with %s\n\n", weightsfile, mapfile);}

if(mode==1||mode==2||mode==3)
{printf("Calculating weights for a second time, using results from first run (stored in %s)\n\n", weightsfile);}
}

///////////////////////////

//checks for genefile/chunks/chunksbp

if(found3>1)
{printf("Error, can only specify one from \"--genefile\", \"--chunks\" and \"--chunks-bp\"\n\n");exit(1);}

if(mode==7&&found3==0)
{printf("Error, to cut into genes it is necessary to provide one from \"--genefile\", \"--chunks\" and \"--chunks-bp\"\n\n");exit(1);}

//won't check if can (will) be used

if(strcmp(genefile,"blank")!=0)	//check genefile exists
{
if((input=fopen(genefile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n",genefile);exit(1);}
fclose(input);
}

///////////////////////////

//check for respfile and sort

if(strcmp(respfile,"blank")!=0)	//check respfile exists
{
if(mode!=9&&mode!=11&&mode!=16&&mode!=17)
{printf("Error, can only provide responses with \"--calc-genes-reml\", \"--calc-genes-bivar\", \"--reml\" or \"--decompose\"\n\n");exit(1);}

if((input=fopen(respfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n", respfile);exit(1);}
fclose(input);

if(mode==9||mode==11||mode==16||mode==17||mode==18)
{
num_resps=countcols(respfile)-2;

if(mpheno2!=-1&&mode!=11&&mode!=18)
{printf("Error, can only use \"--mpheno2\" when performing bivariate analysis\n\n");exit(1);}

if(mode==11||mode==18)
{
if((mpheno!=-1&&mpheno2==-1)||(mpheno==-1||mpheno2!=-1))
{printf("Error, must supply both \"--mpheno\" and \"--mpheno2\"\n\n");exit(1);}
}

if(mpheno!=-1&&mpheno==mpheno2)
{printf("Error, can not provide the same value (%d) for \"--mpheno\" and \"--mpheno2\"\n\n", mpheno);exit(1);}

if(mpheno!=-1&&mpheno>num_resps)
{printf("Error, value %d provided for \"--mpheno\" but %s contains only %d phenotypes\n\n", mpheno, respfile, num_resps);exit(1);}

if(mpheno!=-1&&mpheno2>num_resps)
{printf("Error, value %d provided for \"--mpheno2\" but %s contains only %d phenotypes\n\n", mpheno2, respfile, num_resps);exit(1);}

if((mode==9||mode==16||mode==17)&&num_resps>1&&mpheno==-1)
{printf("Error, %s contains multiple phenotypes, so you must specify one using \"--mpheno\"\n\n", respfile);exit(1);}

if((mode==11||mode==18)&&num_resps>2&&mpheno==-1)
{printf("Error, %s contains more than 2 phenotypes, so you must specify two using \"--mpheno\" and \"--mpheno2\"\n\n", respfile);exit(1);}

if(mpheno==-1&&num_resps==1){mpheno=1;}
if(mpheno==-1&&num_resps==2){mpheno=1;mpheno=2;}

num_resps_use=1;
if(mode==11){num_resps_use=2;}
keepresps[0]=mpheno-1;keepresps[1]=mpheno2-1;

if(num_resps<num_resps_use){printf("Error, %s should contain at least %d columns (not %d)\n\n", respfile, num_resps_use+2, num_resps+2);exit(1);}
}
}

//??? and something here about storing only small pvalues

///////////////////////////

//check for covarfile and sort

if(strcmp(covarfile,"blank")!=0)	//check covarfile exists
{
if(mode!=9&&mode!=11&&mode!=16)
{printf("Error, can only provide covariates with \"--calc-genes-reml\", \"--calc-genes-bivar\" or \"--reml\"\n\n");exit(1);}

if((input=fopen(covarfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n", covarfile);exit(1);}
fclose(input);

num_covars=countcols(covarfile)-1;
}

///////////////////////////

if(mode==19)	//sort out blup
{
if(strcmp(blupfile,"blank")==0&&strcmp(remlfile,"blank")==0)
{printf("Error, to calculate blups must provide either blupfile (and maybe regfile) or remlfile\n\n"); exit(1);}

if(strcmp(regfile,"blank")!=0)	//check regfile exists and get num_regs
{
if((input=fopen(regfile,"r"))==NULL)
{printf("File %s does not exist - correct this first\n\n",regfile);exit(1);}
fclose(input);
if(countcols(regfile)<4)
{printf("Error, %s should have at least 4 columns (not %d)\n\n", regfile, countcols(regfile));exit(1);}
num_regs=countcols(regfile)-4;
}
else
{num_regs=0;}

if(strcmp(blupfile,"blank")!=0)	//check blupfile exists and has correct size
{
if((input=fopen(blupfile,"r"))==NULL)
{printf("File %s does not exist - correct this first\n\n",blupfile);exit(1);}
fclose(input);

if(countcols(blupfile)!=2+2*(num_kins+num_regs))
{printf("Error, if there are %d kinships and %d regions, %s should have %d columns (not %d)\n(To provide kinships, use \"--grm\" or \"--mgrm\"; for regions use \"--regfile\")\n\n", num_kins, num_regs, blupfile, 2+2*(num_kins+num_regs), countcols(blupfile));exit(1);}
}

if(strcmp(remlfile,"blank")!=0)
{
if((input=fopen(remlfile,"r"))==NULL)
{printf("Error, %s does not exist - correct this first\n\n", remlfile);exit(1);}

if(fscanf(input, "Kins:%d Regions:%d ", &readint, &num_regs)!=2)
{printf("Error reading %s; first line should be \"Kins:X Regions:Y\"\n\n",remlfile);exit(1);}

if(readint>0&&num_kins==0)
{printf("Error, necessary to provide kinship stems using either \"--grm\" or \"--mgrm\"\n\n");exit(1);}

if(readint!=num_kins)
{printf("Error, number of kinships provided (%d) does not match number in %s (%d)\n\n", num_kins, remlfile, readint);exit(1);}

if(fscanf(input, "Blupfile: %s ", readstring)!=1)
{printf("Error reading %s; second line should be \"Blupfile: <filename>\"\n\n", remlfile);exit(1);}

if(strcmp(blupfile,"blank")!=0&&strcmp(blupfile,readstring)!=0)
{printf("Error, blupfile provided (%s) does not match that in %s (%s)\n\n", blupfile, remlfile, readstring);exit(1);}
if(strcmp(readstring,"none")!=0){strcpy(blupfile,readstring);}

if(fscanf(input, "Regfile: %s ", readstring)!=1)
{printf("Error reading %s; third line should be \"Regfile: <filename> %s\n\n", remlfile, readstring);exit(1);}

if(strcmp(regfile,"blank")!=0&&strcmp(regfile,readstring)!=0)
{printf("Error, regfile provided (%s) does not match that in %s (%s)\n\n", regfile, remlfile, readstring);exit(1);}
if(strcmp(readstring,"none")!=0){strcpy(regfile,readstring);}

if(strcmp(regfile,"blank")!=0)	//check regfile exists (perhaps again!) and has correct size
{
if((input=fopen(regfile,"r"))==NULL)
{printf("File %s does not exist - correct this first\n\n",regfile);exit(1);}
fclose(input);

if(countcols(regfile)!=4+num_regs)
{printf("Error, %s should have %d columns (not %d). Has file been changed?\n\n", regfile, 4+num_regs, countcols(regfile));exit(1);}
}

if(strcmp(blupfile,"blank")!=0)	//check blupfile exists (perhaps again!) and has correct size
{
if((input=fopen(blupfile,"r"))==NULL)
{printf("File %s does not exist - correct this first\n\n",blupfile);exit(1);}
fclose(input);

if(countcols(blupfile)!=2+2*(num_kins+num_regs))
{printf("Error, %s should have %d columns (not %d). Has file been changed?\n\n", blupfile, 2*(num_kins+num_regs), countcols(blupfile));exit(1);}
}
}

printf("%d kins and %d regs\n", num_kins, num_regs);
}	//end of mode=19

/////////

if(mode==20)
{
if(strcmp(regfile,"blank")==0)
{printf("Error, to calculate scores must provide a scorefile\n\n"); exit(1);}

if(strcmp(regfile,"blank")!=0)	//check scorefile exists and get num_regs
{
if((input=fopen(regfile,"r"))==NULL)
{printf("File %s does not exist - correct this first\n\n",regfile);exit(1);}
fclose(input);
if(countcols(regfile)<4)
{printf("Error, %s should have at least 4 columns (not %d)\n\n", regfile, countcols(regfile));exit(1);}
num_regs=countcols(regfile)-4;
}

printf("%d profiles\n", num_regs);
}

///////////////////////////

//checks for reml priors

if(priora!=-1||priorb!=-1)
{
if(mode!=9)
{printf("Error, can only use \"--priora\" and \"--priorb\" for gene-based reml\n\n");exit(1);}
if(priora==-1||priorb==-1)
{printf("Error, must either provide both \"--priora\" and \"--priorb\" or neither\n\n");exit(1);}
}
else	//so both at -1
{priorb=1.0;}

/////////

if(cut1!=-1||cut2!=-1)
{
if(cut1==-1||cut2==-1)
{printf("Error, to divide into regions, must provide both \"--sig1\" and \"--sig2\" (typical values are 0.000001 and 0.01)\n\n");exit(1);}
}

///////////////////////////

//sort filters

if((mode==25||mode==26||mode==27)&&(minmaf!=-1||maxmaf!=-1||minvar!=-1||minobs!=-1))
{printf("Error, when converting datatypes, can't use any of \"--minmaf\", \"--maxmaf\", \"--minvar\" or \"--minobs\"\n\n");exit(1);}

if(threshold!=-1&&((mode!=25&&mode!=26&&mode!=27)||strcmp(chiamofile, "blank")!=0||mode!=26))
{printf("Error, can only use \"--threshold\" when converting datatypes from chiamo to binary PLINK\n\n");exit(1);}

if(minmaf==-1)
{
if(mode==25||mode==26||mode==27){minmaf=0;}
else{minmaf=0.01;}
}
if(maxmaf==-1)
{
maxmaf=0.5;
}
if(minvar==-1)
{
if(mode==25||mode==26||mode==27){minvar=0;}
else{minvar=0.01;}
}
if(minobs==-1)
{
if(mode==25||mode==26||mode==27){minobs=0;}
else{minobs=0.95;}
}

if(threshold==-1)
{
if(mode==25||mode==26||mode==27){threshold=0.95;}
}



