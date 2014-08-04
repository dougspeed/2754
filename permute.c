/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


int permute_int(int *array, int length)
//permutes array of size length
{
int i, j;
int temp;

for(i=1;i<length;i++)
{
j=rand()%length;
temp=array[i];
array[i]=array[j];
array[j]=temp;
}

return(0);
}	//end of permute



int permute_double(double *array, int length)
//permutes array of size length
{
int i, j;
double temp;

for(i=1;i<length;i++)
{
j=rand()%length;
temp=array[i];
array[i]=array[j];
array[j]=temp;
}

return(0);
}	//end of permute_double







