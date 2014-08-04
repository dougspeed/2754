/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


//first set the filename of the details file
if(mode==2||mode==3)
{
if(strcmp(weightsfile,"blank")==0){sprintf(filename,"%ssection_details.txt", folder);}
else{sprintf(filename,"%sre-section_details.txt", folder);}
}
if(mode==5||mode==6){sprintf(filename,"%spartition_details.txt", folder);}
if(mode==8||mode==9||mode==10||mode==11||mode==12){sprintf(filename,"%sgene_details.txt", folder);}

//and open
if((input=fopen(filename,"r"))==NULL)
{
if(mode==2||mode==3){printf("Can't find file %s which contains details of sections\nThis will have been created using argument \"--cut-weights\"\n\n", filename);exit(1);}
if(mode==5||mode==6){printf("Can't find file %s which contains details of regions\nThis will have been created using argument \"--cut-kins\"\n\n", filename);exit(1);}
if(mode==8||mode==9||mode==10||mode==11||mode==12){printf("Can't find file %s which contains details of genes\nThis will have been created using argument \"--cut-genes\"\n\n", filename);exit(1);}
}

//now check (get) mapfile from the top line
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading top line of %s. Has file been changed since creation with argument \"--cut-weights\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafiles:")!=0)
{printf("Error reading top line of %s. Has file been changed since creation with argument \"--cut-weights\"\n\n", filename);exit(1);}

if(strcmp(mapfile,"blank")!=0&&strcmp(mapfile,readstring3)!=0)
{
if(mode==2||mode==3)
{printf("Mapfile used with argument \"--cut-weights\" (%s) is different from that supplied now (%s)\n\n", readstring3, mapfile);exit(1);}
if(mode==5||mode==6)
{printf("Mapfile used with argument \"--cut-kins\" (%s) is different from that supplied now (%s)\n\n", readstring3, mapfile);exit(1);}
if(mode==8||mode==9||mode==10||mode==11||mode==12)
{printf("Mapfile used with argument \"--cut-genes\" (%s) is different from that supplied now (%s)\n\n", readstring3, mapfile);exit(1);}
}
else
{
if(mode==3||mode==6||mode==10||mode==11||mode==12){strcpy(mapfile,readstring3);}
}

//now read second line to see if using predictor keepfile
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading second line of %s. Has file been changed since creation with argument \"--cut-weights\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Using")!=0&&strcmp(readstring,"Extracting")!=0&&strcmp(readstring,"Calculating")!=0)
{printf("Error reading second line of %s. Has file been changed since creation with argument \"--cut-weights\"\n\n", filename);exit(1);}

if(strcmp(bpredfile,"blank")!=0)	//keepfile provided, expect to see same one
{
if(strcmp(readstring,"Extracting")!=0)	//so none provided when cutting
{
if(mode==2||mode==3)
{printf("A predictor file is provided now (%s), but this was not provided when cutting with argument \"--cut-weights\n\n", bpredfile);exit(1);}
if(mode==5||mode==6)
{printf("A predictor file is provided now (%s), but this was not provided when cutting with argument \"--cut-kins\n\n", bpredfile);exit(1);}
if(mode==8||mode==9||mode==10||mode==11||mode==12)
{printf("A predictor file is provided now (%s), but this was not provided when cutting with argument \"--cut-genes\n\n", bpredfile);exit(1);}
}
if(strcmp(readstring3,bpredfile)!=0)	//so different one provided when cutting
{
if(mode==2||mode==3)
{printf("A predictor file is provided now (%s), but a different one was provided when cutting with argument \"--cut-weights (%s)\n\n", bpredfile, readstring3);exit(1);}
if(mode==5||mode==6)
{printf("A predictor file is provided now (%s), but a different one was provided when cutting with argument \"--cut-kins (%s)\n\n", bpredfile, readstring3);exit(1);}
if(mode==8||mode==9||mode==10||mode==11||mode==12)
{printf("A predictor file is provided now (%s), but a different one was provided when cutting with argument \"--cut-genes (%s)\n\n", bpredfile, readstring3);exit(1);}
}
}
else		//so bpredfile not used
{
if(strcmp(readstring,"Extracting")==0)		//but one provided when cutting
{
if(mode==2)
{printf("The predictor file %s was specified when using argument \"--cut-weights\"\n\n", readstring3);exit(1);}
if(mode==5)
{printf("The predictor file %s was specified when using argument \"--cut-kins\"\n\n", readstring3);exit(1);}
if(mode==8||mode==9||mode==15)
{printf("The predictor file %s was specified when using argument \"--cut-genes\"\n\n", readstring3);exit(1);}
//only here if model 3, 6, 10 or 16
strcpy(bpredfile,readstring3);
}
}

fclose(input);











