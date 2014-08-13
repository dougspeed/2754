/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


data=malloc(sizeof(double)*num_samples_use*bitsize);
if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0)
{datainput=open_data_fly(chiamofile, spfile, num_samples, chiamoprobs, chiamoheaders);}

//set filenames for data and open
if(mode==25)	//sp
{
sprintf(filename,"%s_out.sp", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s\n\n",filename);exit(1);}
}
if(mode==26)	//bed
{
sprintf(filename,"%s_out.bed", outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename);exit(1);}
gen=108;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=27;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=1;fwrite(&gen, sizeof(unsigned char), 1, output);
}
if(mode==27)	//speed
{
sprintf(filename,"%s_out.speed", outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s\n\n",filename);exit(1);}
}

bittotal=(int)((data_length-1)/bitsize+1);
current=0;weightsum=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

current=read_data_fly(bedfile, chiamofile, spfile, speedfile, datainput, current, data, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, num_samples, num_preds, missingvalue, chiamoheaders, chiamoprobs, al1, al2, -1);

for(j=0;j<bitlength;j++)
{
gen=0;
for(i=0;i<num_samples_use;i++)
{
value=data[i+j*num_samples_use];

if(mode==25)
{
if(value!=missingvalue)
{
if(value==(int)value){fprintf(output,"%d ", (int)value);}
else{fprintf(output,"%f ", value);}
}
else
{fprintf(output,"NA ");}
}
if(mode==26)
{
gen2=1;
if(value<=0.1){gen2=3;}
if(value>=1.9){gen2=0;}
if(value>0.95&&value<1.05){gen2=2;}
if(j==0&&i<10){printf("%d %d\n", i+1, gen2);}
gen+=gen2*pow(2,2*(i%4));
if(i%4==3||i==num_samples_use-1){fwrite(&gen, sizeof(unsigned char), 1, output);gen=0;}
}
if(mode==27)
{
fwrite(&value, sizeof(float), 1, output);
}
}	//end of i loop
if(mode==25){fprintf(output,"\n");}
}	//end of j loop
}	//end of bit loop

if(strcmp(chiamofile,"blank")!=0||strcmp(spfile,"blank")!=0){fclose(datainput);}
free(data);

//write bimfile
sprintf(filename2,"%s_out.bim", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s\n\n",filename2);exit(1);}
for(j=0;j<data_length;j++)
{fprintf(output2, "%d\t%s\t0\t%ld\t%c\t%c\n", chr[j], prednames[j], (long int)bp[j], al1[j], al2[j]);}
fclose(output2);

//write famfile
sprintf(filename3,"%s_out.fam", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s\n\n",filename3);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output3, "%s %s 0 0 0 0\n", ids1[i], ids2[i]);}
fclose(output3);


