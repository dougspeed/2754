/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


if(mode==2&&section==-1)
{printf("Necessary to provide section number to consider using argument \"--section\"\n\n");exit(1);}

if(mode==1||mode==2||mode==3)
{printf("If you plan to perform \"Genomic Partitioning\", make sure you calculate weightings over the correct predictors (the union of all regions.\nSee http://dougspeed.com/ldak/genomic-partitioning for details.\n\n");}

if(mode==5||mode==8||mode==9||mode==11)
{
if(partition==-1){printf("Error, necessary to provide partition number to consider using argument \"--partition\"\n\n");exit(1);}
}

if(mode==5||mode==7||mode==8||mode==9||mode==11||(mode==16&&num_regs!=0)||(mode==24&&num_regs>0))
{
if(strcmp(weightsfile,"blank")==0&&ignoreweights==0){printf("Necessary to provide weightsfile using argument \"--weights\" or for uniform weights use argument \"--ignore-weights YES\"\n\n");exit(1);}
}

if(mode==9||mode==11||mode==16)
{
if(strcmp(respfile,"blank")==0)
{printf("Error, necessary to provide phenotypic values using argument \"--pheno\"\n\n");exit(1);}
}

if(mode==17&&num_kins!=1)
{printf("Error, necessary to provide a (single) kinship file using either \"--grm\" or \"--mgrm\"\n\n");exit(1);}

if(mode==21&&her==-1)
{printf("Error, necessary to provide heritability using argument \"--her\"\n\n");exit(1);}

if(mode==22&&num_samples==-1)
{printf("Error, necessary to provide number of samples using argument \"--num-samples\"\n\n");exit(1);}

if(mode==22&&num_snps==-1)
{printf("Error, necessary to provide number of SNPs using argument \"--num-snps\"\n\n");exit(1);}

if(mode==23||mode==24)
{
if(strcmp(kinlist,"blank")==0&&strcmp(kinname,"blank")==0)
{printf("Error, necessary to provide kinship stems using either \"--grm\" or \"--mgrm\"\n\n");exit(1);}

if(strcmp(bsampfile,"blank")!=0)
{printf("Error, can not keep individuals when adding or subtracting kinships\n\n");exit(1);}

if(mode==24&&num_kins==1&&num_regs==0)
{printf("Warning, with argument \"--sub-grm\" it is usual to supplie either more than one kinship files or regions to subtract\n\n");}
}


