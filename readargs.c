/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int print1()
{
printf("Can only specify one argument from \"--cut-weights\", \"--calc-weights\", \"--join-weights\", \"--cut-kins\",\"--calc-kins\", \"--join-kins\", \"--cut-genes\", \"--calc-genes-kin\", \"--calc-genes-reml\", \"--join-genes-reml\", \"--reml\", \"--decompose\", \"--calc-blups\", \"--make-phenos\", \"--make-snps\", \"--add-grm\", \"--sub-grm\", \"--make-sp\", \"--make-bed\" or \"--make-speed\"\n\n");
return(0);
}


count=1;found2=0;found3=0;
while(count<argc)
{
strcpy(argname,argv[count]);

if(argc==count+1){printf("Argument \"%s\" does not have corresponding value - arguments must be in pairs\n\n",argname);exit(1);}

//see which argument specified, and get the corresponding number/file
found=0;

/////////

if(strcmp(argname,"--cut-weights")==0)
{
if(mode!=0){print1();exit(1);}
mode=1;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-weights")==0)
{
if(mode!=0){print1();exit(1);}
mode=2;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--join-weights")==0)
{
if(mode!=0){print1();exit(1);}
mode=3;strcpy(folder3,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--cut-kins")==0)
{
if(mode!=0){print1();exit(1);}
mode=4;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-kins")==0)
{
if(mode!=0){print1();exit(1);}
mode=5;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--join-kins")==0)
{
if(mode!=0){print1();exit(1);}
mode=6;strcpy(folder3,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--cut-genes")==0)
{
if(mode!=0){print1();exit(1);}
mode=7;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-genes-kins")==0)
{
if(mode!=0){print1();exit(1);}
mode=8;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-genes-reml")==0)
{
if(mode!=0){print1();exit(1);}
mode=9;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--join-genes-reml")==0)
{
if(mode!=0){print1();exit(1);}
mode=10;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-genes-bivar")==0)
{
if(mode!=0){print1();exit(1);}
mode=11;strcpy(folder3,argv[count+1]);found=1;}

if(strcmp(argname,"--join-genes-bivar")==0)
{
if(mode!=0){print1();exit(1);}
mode=12;strcpy(folder3,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--reml")==0)
{
if(mode!=0){print1();exit(1);}
mode=16;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--decompose")==0)
{
if(mode!=0){print1();exit(1);}
mode=17;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-blups")==0)
{
if(mode!=0){print1();exit(1);}
mode=19;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--calc-scores")==0)
{
if(mode!=0){print1();exit(1);}
mode=20;strcpy(outfile2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--make-phenos")==0)
{
if(mode!=0){print1();exit(1);}
mode=21;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--make-snps")==0)
{
if(mode!=0){print1();exit(1);}
mode=22;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--add-grm")==0)
{
if(mode!=0){print1();exit(1);}
mode=23;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--sub-grm")==0)
{
if(mode!=0){print1();exit(1);}
mode=24;strcpy(outfile2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--make-sp")==0)
{
if(mode!=0){print1();exit(1);}
mode=25;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--make-bed")==0)
{
if(mode!=0){print1();exit(1);}
mode=26;strcpy(outfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--make-speed")==0)
{
if(mode!=0){print1();exit(1);}
mode=27;strcpy(outfile2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--workdir")==0)
{strcpy(workdir2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--bfile")==0)
{
sprintf(bedfile2,"%s.bed", argv[count+1]);
sprintf(mapfile2,"%s.bim", argv[count+1]);
sprintf(famfile2,"%s.fam", argv[count+1]);
sprintf(datafile2,"%s.bed",argv[count+1]);
found=1;found2++;}

if(strcmp(argname,"--chiamo")==0)
{
sprintf(chiamofile3,"%s", argv[count+1]);
sprintf(mapfile2,"%s.map", argv[count+1]);
sprintf(famfile2,"%s.fam", argv[count+1]);
sprintf(datafile3,"%s",argv[count+1]);
found=1;found2++;}

if(strcmp(argname,"--sp")==0)
{
sprintf(spfile2,"%s.sp", argv[count+1]);
sprintf(mapfile2,"%s.bim", argv[count+1]);
sprintf(famfile2,"%s.fam", argv[count+1]);
sprintf(datafile2,"%s.sp",argv[count+1]);
found=1;found2++;}

if(strcmp(argname,"--speed")==0)
{
sprintf(speedfile2,"%s.speed", argv[count+1]);
sprintf(mapfile2,"%s.bim", argv[count+1]);
sprintf(famfile2,"%s.fam", argv[count+1]);
sprintf(datafile2,"%s.speed",argv[count+1]);
found=1;found2++;}

/////////

if(strcmp(argname,"--chiamo-probs")==0)
{chiamoprobs=atoi(argv[count+1]);found=1;
if(chiamoprobs<2||chiamoprobs>4){printf("%s should be followed by 2, 3, or 4 (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--chiamo-headers")==0)
{chiamoheaders=atoi(argv[count+1]);found=1;
if(chiamoheaders<2){printf("%s should be followed by an integer >=2 (not %s)\n(chimaofile headers must (at least) provide two alleles)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--chiamo-suffix")==0)
{strcpy(chiamosuffix,argv[count+1]);found=1;}

if(strcmp(argname,"--missing-value")==0)
{missingvalue=atof(argv[count+1]);found=1;}

if(strcmp(argname,"--power")==0)
{power=atof(argv[count+1]);found=1;
if(power<-2||power>1){printf("%s can be followed by any float, but unusual to be outside the range [-2,1] (currently %s), default -1\n\n", argname, argv[count+1]);exit(1);}
}

/////////

if(strcmp(argname,"--extract-index")==0)
{printf("Sorry, argument \"--extract-index\" no longer allowed; provide predictor names using \"--extract\" instead\n\n");}

if(strcmp(argname,"--extract")==0)
{strcpy(bpredfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--keep-index")==0)
{printf("Sorry, argument \"--keep-index\" no longer allowed; provide sample ids using \"--keep\" instead\n\n");}

if(strcmp(argname,"--keep")==0)
{strcpy(bsampfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--subset-number")==0)
{num_subs=atoi(argv[count+1]);found=1;
if(num_subs==1){printf("%s should not be followed by 1\n(If determined to use only one subset, use \"--keep-index\" instead)\n\n", argname);exit(1);}
if(num_subs<=1){printf("%s should be followed by a positive integer >1 (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--subset-prefix")==0)
{strcpy(subpref2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--minmaf")==0)
{minmaf=atof(argv[count+1]);found=1;
if(minmaf<0||minmaf>=0.5){printf("%s should be followed by a float in [0,0.5) (not %s), default 0.01\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--maxmaf")==0)
{maxmaf=atof(argv[count+1]);found=1;
if(maxmaf<=0||maxmaf>0.5){printf("%s should be followed by a float in (0,.5] (not %s), default 0.01\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--minvar")==0)
{minvar=atof(argv[count+1]);found=1;
if(minvar<=0){printf("%s should be followed by a positive float (not %s), default 0.01\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--minobs")==0)
{minobs=atof(argv[count+1]);found=1;
if(minobs<=0||minobs>1){printf("%s should be followed by a float in (0,1] (not %s), default 0.01\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--encoding")==0)
{encoding=atoi(argv[count+1]);found=1;
if(encoding!=1&&encoding!=2&&encoding!=3&&encoding!=4){printf("%s should be followed by either 1 (additive), 2 (dominant), 3 (recessive) or 4 (heterozygotes) (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--threshold")==0)
{threshold=atof(argv[count+1]);found=1;
if(threshold<=0.5||threshold>1){printf("%s should be followed by a float in (0.5,1] (not %s), default 0.95\n\n", argname, argv[count+1]);exit(1);}
}

///////////////////////////

if(strcmp(argname,"--weights")==0)
{strcpy(weightsfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--ignore-weights")==0)
{
ignoreweights=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){ignoreweights=1;}
if(strcmp(argv[count+1],"NO")==0){ignoreweights=0;}
if(ignoreweights==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

/////////

if(strcmp(argname,"--section")==0)
{section=atoi(argv[count+1]);found=1;
if(section<=0){printf("%s should be followed by a positive integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--section-length")==0)
{section_length=atoi(argv[count+1]);found=1;
if(section_length<=0){printf("%s should be followed by a positive integer (not %s), default 3000 predictors\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--section-buffer")==0)
{section_buffer=atoi(argv[count+1]);found=1;
if(section_buffer<0){printf("%s should be followed by a non-negative integer (not %s), default 500 predictors\n\n", argname, argv[count+1]);exit(1);}
}

/////////

if(strcmp(argname,"--partition")==0)
{partition=atoi(argv[count+1]);found=1;
if(partition<=0){printf("%s should be followed by a positive integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--partition-length")==0)
{partition_length=atoi(argv[count+1]);found=1;
if(partition_length<=0){printf("%s should be followed by a positive integer (not %s), default 500000\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--partition-number")==0)
{num_parts=atoi(argv[count+1]);found=1;
if(num_parts<0){printf("%s should be followed by a non-negative integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--partition-prefix")==0)
{strcpy(partpref2,argv[count+1]);found=1;}

if(strcmp(argname,"--by-chr")==0)
{
bychr=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){bychr=1;}
if(strcmp(argv[count+1],"NO")==0){bychr=0;}
if(bychr==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--kinship-matrix")==0)
{
kinraw=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kinraw=1;}
if(strcmp(argv[count+1],"NO")==0){kinraw=0;}
if(kinraw==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--kinship-gz")==0)
{
kingz=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kingz=1;}
if(strcmp(argv[count+1],"NO")==0){kingz=0;}
if(kingz==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

/////////

if(strcmp(argname,"--region-number")==0)
{num_regs=atoi(argv[count+1]);found=1;
if(num_regs<0){printf("%s should be followed by a non-negative integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--region-prefix")==0)
{strcpy(regpref2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--genefile")==0)
{strcpy(genefile2,argv[count+1]);found=1;found3++;}

if(strcmp(argname,"--chunks")==0)
{chunks=atof(argv[count+1]);found=1;found3++;
if(chunks<=0){printf("%s should be followed by +tive float (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--chunks-bp")==0)
{chunksbp=atoi(argv[count+1]);found=1;found3++;
if(chunksbp<=0){printf("%s should be followed by +tive integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--gene-buffer")==0)
{gene_buffer=atoi(argv[count+1]);found=1;
if(gene_buffer<0){printf("%s should be non-negative integer (not %s), default 0 predictors\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--min-weight")==0)
{minweight=atof(argv[count+1]);found=1;
if(minweight<0){printf("%s should be followed by +tive float or zero (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--overlap")==0)
{overlap=-1;found=1;
if(strcmp(argv[count+1],"NO")==0){overlap=0;}
if(strcmp(argv[count+1],"YES")==0){overlap=1;}
if(overlap==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--pvafile")==0)
{strcpy(pvafile2,argv[count+1]);found=1;}

////////

if(strcmp(argname,"--pheno")==0)
{strcpy(respfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--mpheno")==0)
{mpheno=atoi(argv[count+1]);found=1;
if(mpheno<0){printf("%s should be non-negative integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--mpheno2")==0)
{mpheno2=atoi(argv[count+1]);found=1;
if(mpheno2<0){printf("%s should be non-negative integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--permute")==0)
{
permute=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){permute=1;}
if(strcmp(argv[count+1],"NO")==0){permute=0;}
if(permute==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--covar")==0)
{strcpy(covarfile2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--grm")==0)
{strcpy(kinname2,argv[count+1]);found=1;}

if(strcmp(argname,"--mgrm")==0)
{strcpy(kinlist2,argv[count+1]);found=1;}

if(strcmp(argname,"--kinship-details")==0)
{
kindetails=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){kindetails=1;}
if(strcmp(argv[count+1],"NO")==0){kindetails=0;}
if(kindetails==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--eigen")==0)
{strcpy(eigenfile2,argv[count+1]);found=1;}

/////////

if(strcmp(argname,"--remlfile")==0)
{strcpy(remlfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--blupfile")==0)
{strcpy(blupfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--regfile")==0)
{strcpy(regfile2,argv[count+1]);found=1;}

if(strcmp(argname,"--scorefile")==0)
{strcpy(regfile2,argv[count+1]);found=1;}

////////

if(strcmp(argname,"--num-phenos")==0)
{num_phenos=atoi(argv[count+1]);found=1;
if(num_phenos<=0){printf("%s should be followed by +tive integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--num-causals")==0)
{num_causals=atoi(argv[count+1]);found=1;
if(num_causals<=0&&num_causals!=-1){printf("%s should be followed by +tive integer (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--her")==0)
{her=atof(argv[count+1]);found=1;
if(her<=0||her>1){printf("%s should be followed by float within (0,1] (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--causals")==0)
{strcpy(causalsfile2,argv[count+1]);found=1;
}

if(strcmp(argname,"--effects")==0)
{strcpy(effectsfile2,argv[count+1]);found=1;
}

/////////

if(strcmp(argname,"--priora")==0)
{priora=atof(argv[count+1]);found=1;
if(priora<=0){printf("%s should be followed by positive float (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--priorb")==0)
{priorb=atof(argv[count+1]);found=1;
if(priorb<=0){printf("%s should be followed by positive float (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--gwide-adjust")==0)
{
adjust=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){adjust=1;}
if(strcmp(argv[count+1],"NO")==0){adjust=0;}
if(adjust==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--sig1")==0)
{cut1=atof(argv[count+1]);found=1;
if(cut1<=0){printf("%s should be followed by float within (0,1) (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--sig2")==0)
{cut2=atof(argv[count+1]);found=1;
if(cut2<0){printf("%s should be followed by float within [0,1) (not %s)\n\n", argname, argv[count+1]);exit(1);}
}

if(strcmp(argname,"--score-test")==0)
{
scoretest=-1;found=1;
if(strcmp(argv[count+1],"YES")==0){scoretest=1;}
if(strcmp(argv[count+1],"NO")==0){scoretest=0;}
if(scoretest==-1){printf("%s should be followed by YES or NO (not %s)\n\n", argname, argv[count+1]);exit(1);}
}



if(found==0){printf("%s is not a recognised argument\n",argname);exit(1);}
count+=2;
}


printf("Arguments:\n");
count=1;
while(count<argc)
{printf("%s %s\n", argv[count], argv[count+1]);count+=2;}
printf("\n");


if(mode==0)
{
printf("At least one pair of arguments is required, including ONE main argument\n\n");

printf("Main argument must be one of:\n\n");
printf("--cut-weights <folder> (requires bfile/chiamo/sp/speed) - break genome into sections for calculating weights\n\
--calc-weights <folder> (requires bfile/chiamo/sp/speed and section) - calculate weights for section specified\n\
--join-weights <folder> - join up weights\n\n");

printf("--cut-kins <folder> (requires bfile/chiamo/sp/speed) - break genome into partitions for calculating kinships\n\
--calc-kins <folder> (requires bfile/chiamo/sp/speed, region and weights or ignoreweights YES) - calculate kinship for partition specified\n\
--join-kins <folder> - join up kinships\n\n");

printf("--cut-genes <folder> (requires bfile/chiamo/sp/speed and genefile, chunks or chunks-bp) - calculate breakpoints and break genome into groups\n\
--calc-genes-kins <folder> (requires bfile/chiamo/sp/speed, group and weights or ignoreweights YES) - calculate kinships for each gene in groups\n\
--calc-genes-reml <folder> (requires bfile/chiamo/sp/speed, group and weights or ignoreweights YES) - performs regression for each gene in groups\n\
--join-genes-reml <folder> - joins together regression results for all groups\n\n");

printf("--reml <output> (requires grm, mgrm or region-number/region-prefix and pheno) - performs generalised reml analysis\n\
--decompose <output> (requires grm) - computes the eigen-decomposition\n\n");

printf("--make-phenos <outfile> (requires bfile/chiamo/sp/speed and her) - simulate phenotypes\n\
--make-snps <outfile> (requires num-samples and num-snps) - simulate SNPs\n\n");

printf("--add-grm <outfile> (requires mgrm) - add kinships\n\
--sub-grm <outfile> (requires mgrm) - subtract kinships\n\n");

printf("--make-sp <outfile> (requires bfile/chiamo/speed) - converts to sp format\n\
--make-bed <outfile> (requires chiamo/sp/speed) - converts to bed format\n\
--make-speed <outfile> (requires bfile/chiamo/sp) - converts to sp binary format\n\n");

exit(1);
}




