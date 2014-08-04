/*
Copyright 2014 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.r

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/


/* The ziggurat method for RNOR and REXP
Combine the code below with the main program in which you want
normal or exponential variates.   Then use of RNOR in any expression
will provide a standard normal variate with mean zero, variance 1,
while use of REXP in any expression will provide an exponential variate
with density exp(-x),x>0.
Before using RNOR or REXP in your main, insert a command such as
zigset(86947731 );
with your own choice of seed value>0, rather than 86947731.
(If you do not invoke zigset(...) you will get all zeros for RNOR and REXP.)
For details of the method, see Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software.
*/

//HAVE MADE ONE EDIT IN LINE 36 (long to int) TO MAKE TRANSFERABLE TO 64-BIT MACHINES (DOUG SPEED)

/* Also have cdf function for standard normal
Code obtained from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
based on Abromowitz and Stegun approximation approximation explained
http://www.math.sfu.ca/~cbm/aands/page_932.htm
*/


static unsigned int jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static int hz;
static unsigned int iz, kn[128], ke[256];
static float wn[128],fn[128],we[256],fe[256];

#define RNOR (hz=SHR3, iz=hz&127, (fabs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, (    jz <ke[iz])? jz*we[iz] : efix())

/* nfix() generates variates from the residue when rejection in RNOR occurs. */

float nfix(void)
{
const float r = 3.442620f;     /* The start of the right tail */
static float x, y;
 for(;;)
  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }
                         /* iz>0, handle the wedges of other strips */
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

     /* initiate, try to exit for(;;) for loop*/
      hz=SHR3;
      iz=hz&127;
      if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
  }

}

/* efix() generates variates from the residue when rejection in REXP occurs. */
float efix(void)
{ float x;
 for(;;)
  {  if(iz==0) return (7.69711-log(UNI));          /* iz==0 */
     x=jz*we[iz]; if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);

      /* initiate, try to exit for(;;) loop */
   jz=SHR3;
   iz=(jz&255);
   if(jz<ke[iz]) return (jz*we[iz]);
  }
}
/*--------This procedure sets the seed and creates the tables------*/

void zigset(unsigned long jsrseed)
{  const double m1 = 2147483648.0, m2 = 4294967296.;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
   int i;
   jsr^=jsrseed;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }

/* Set up tables for REXP */
    q = ve/exp(-de);
    ke[0]=(de/q)*m2;
    ke[1]=0;

    we[0]=q/m2;
    we[255]=de/m2;

    fe[0]=1.;
    fe[255]=exp(-de);

   for(i=254;i>=1;i--)
  {de=-log(ve/de+exp(-de));
   ke[i+1]= (de/te)*m2;
   te=de;
   fe[i]=exp(-de);
   we[i]=de/m2;
  }
}



////////////////////////////////////////
//cdf of Normal code

double cdfN(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}



