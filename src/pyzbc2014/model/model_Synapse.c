/*
  This file is part of a variant of version 5.2 of the Zilany, Bruce, and Carney (2014)
  auditory-nerve model. Licensing information can be found below. See README.md for a list
  of changes made versus the original version, which is available at:
  https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx.


  Papers related to this code are cited below:

  Zilany, M. S., & Bruce, I. C. (2006). Modeling auditory-nerve responses for high sound 
  pressure levels in the normal and impaired auditory periphery. The Journal of the 
  Acoustical Society of America, 120(3), 1446-1466.

  * Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A
  Phenomenological model of the synapse between the inner hair cell and auditory
  nerve : Long-term adaptation with power-law dynamics," Journal of the
  Acoustical Society of America 126(5): 2390-2412.
 
  Ibrahim, R. A., and Bruce, I. C. (2010). "Effects of peripheral tuning
  on the auditory nerve's representation of speech envelope and temporal fine
  structure cues," in The Neurophysiological Bases of Auditory Perception, eds.
  E. A. Lopez-Poveda and A. R. Palmer and R. Meddis, Springer, NY, pp. 429-438.

  Zilany, M.S.A., Bruce, I.C., Ibrahim, R.A., and Carney, L.H. (2013).
  "Improved parameters and expanded simulation options for a model of the
  auditory periphery," in Abstracts of the 36th ARO Midwinter Research Meeting.

  * Zilany, M. S., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and expanded 
  simulation options for a model of the auditory periphery. The Journal of the Acoustical 
  Society of America, 135(1), 283-286.

  Please cite these papers marked with asterisks (*) if you publish any research results 
  obtained with this code or any modified versions of this code.


  Zilany, Bruce, Carney (2014) auditory-nerve model source code, Copyright (C) 2024 
    M. S. A. Zilany <msazilany@gmail.com> 
    Ian C. Bruce <ibruce@ieee.org> 
    Laurel H. Carney <laurel_carney@urmc.rochester.edu>
    Daniel R. Guest <daniel_guest@urmc.rochester.edu>

  This program is free software: you can redistribute it and/or modify it under the terms of
  the GNU Affero General Public License as published by the Free Software Foundation, either
  version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License along with this
  program. If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "complex.h"

#define MAXSPIKES 1000000
#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

double Synapse(
    double *ihcout, 
    double *randNums, 
    double tdres, 
    double cf, 
    int totalstim, 
    int nrep,
    double spont, 
    double implnt, 
    double sampFreq,
    double *synouttmp 
) {
    /* Initalize Variables */
    int z, b;
    int resamp = (int)ceil(1 / (tdres * sampFreq));
    double incr = 0.0;
    int delaypoint = (int)floor(7500 / (cf / 1e3));

    double alpha1, beta1, I1, alpha2, beta2, I2, binwidth;
    int k, j, indx, i;
    double synstrength, synslope, CI, CL, PG, CG, VL, PL, VI;
    double cf_factor, PImax, kslope, Ass, Asp, TauR, TauST, Ar_Ast, PTS, Aon,
        AR, AST, Prest, gamma1, gamma2, k1, k2;
    double VI0, VI1, alpha, beta, theta1, theta2, theta3, vsat, tmpst, tmp, PPI,
        CIlast, temp;

    double *sout1, *sout2, *synSampOut, *powerLawIn, *exponOut, *TmpSyn;
    double *m1, *m2, *m3, *m4, *m5;
    double *n1, *n2, *n3;

    double *sampihc, *ihcdims;

    exponOut = (double*)calloc((long) ceil(totalstim*nrep),sizeof(double));
    powerLawIn = (double*)calloc((long) ceil(totalstim*nrep+3*delaypoint),sizeof(double));
    sout1 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    sout2 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    synSampOut  = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    TmpSyn  = (double*)calloc((long) ceil(totalstim*nrep+2*delaypoint),sizeof(double));

    m1 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m2 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m3  = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m4 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m5  = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));

    n1 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    n2 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    n3 = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));

    /*----------------------------------------------------------*/
    /*------- Parameters of the Power-law function -------------*/
    /*----------------------------------------------------------*/
    binwidth = 1/sampFreq;
    alpha1 = 2.5e-6*100e3; beta1 = 5e-4; I1 = 0;
    alpha2 = 1e-2*100e3; beta2 = 1e-1; I2 = 0;

    /*----------------------------------------------------------*/
    /*----- Double Exponential Adaptation ----------------------*/
    /*----------------------------------------------------------*/
    if (spont==100) cf_factor = __min(800,pow(10,0.29*cf/1e3 + 0.7));
    if (spont==4)   cf_factor = __min(50,2.5e-4*cf*4+0.2);
    if (spont==0.1) cf_factor = __min(1.0,2.5e-4*cf*0.1+0.15);

    PImax  = 0.6;                /* PI2 : Maximum of the PI(PI at steady state) */
    kslope = (1+50.0)/(5+50.0)*cf_factor*20.0*PImax;
    /* Ass    = 300*TWOPI/2*(1+cf/100e3); */  /* Older value: Steady State Firing Rate eq.10 */
    Ass    = 800*(1+cf/100e3);    /* Steady State Firing Rate eq.10 */

    if (implnt==1) Asp = spont*3.0;   /* Spontaneous Firing Rate if actual implementation */
    if (implnt==0) Asp = spont*2.75; /* Spontaneous Firing Rate if approximate implementation */
    TauR   = 2e-3;               /* Rapid Time Constant eq.10 */
    TauST  = 60e-3;              /* Short Time Constant eq.10 */
    Ar_Ast = 6;                  /* Ratio of Ar/Ast */
    PTS    = 3;                  /* Peak to Steady State Ratio, characteristic of PSTH */

    /* now get the other parameters */
    Aon    = PTS*Ass;                          /* Onset rate = Ass+Ar+Ast eq.10 */
    AR     = (Aon-Ass)*Ar_Ast/(1+Ar_Ast);      /* Rapid component magnitude: eq.10 */
    AST    = Aon-Ass-AR;                       /* Short time component: eq.10 */
    Prest  = PImax/Aon*Asp;                    /* eq.A15 */
    CG  = (Asp*(Aon-Asp))/(Aon*Prest*(1-Asp/Ass));    /* eq.A16 */
    gamma1 = CG/Asp;                           /* eq.A19 */
    gamma2 = CG/Ass;                           /* eq.A20 */
    k1     = -1/TauR;                          /* eq.8 & eq.10 */
    k2     = -1/TauST;                         /* eq.8 & eq.10 */
            /* eq.A21 & eq.A22 */
    VI0    = (1-PImax/Prest)/(gamma1*(AR*(k1-k2)/CG/PImax+k2/Prest/gamma1-k2/PImax/gamma2));
    VI1    = (1-PImax/Prest)/(gamma1*(AST*(k2-k1)/CG/PImax+k1/Prest/gamma1-k1/PImax/gamma2));
    VI  = (VI0+VI1)/2;
    alpha  = gamma2/k1/k2;       /* eq.A23,eq.A24 or eq.7 */
    beta   = -(k1+k2)*alpha;     /* eq.A23 or eq.7 */
    theta1 = alpha*PImax/VI;
    theta2 = VI/PImax;
    theta3 = gamma2-1/PImax;

    PL  = ((beta-theta2*theta3)/theta1-1)*PImax;  /* eq.4' */
    PG  = 1/(theta3-1/PL);                        /* eq.5' */
    VL  = theta1*PL*PG;                           /* eq.3' */
    CI  = Asp/Prest;                              /* CI at rest, from eq.A3,eq.A12 */
    CL  = CI*(Prest+PL)/PL;                       /* CL at rest, from eq.1 */

    if(kslope>=0)  vsat = kslope+Prest;
    tmpst  = log(2)*vsat/Prest;
    if(tmpst<400) synstrength = log(exp(tmpst)-1);
    else synstrength = tmpst;
    synslope = Prest/log(2)*synstrength;

    k = 0;
    for (indx=0; indx<totalstim*nrep; ++indx)
    {
        tmp = synstrength*(ihcout[indx]);
        if(tmp<400) tmp = log(1+exp(tmp));
        PPI = synslope/synstrength*tmp;

        CIlast = CI;
        CI = CI + (tdres/VI)*(-PPI*CI + PL*(CL-CI));
        CL = CL + (tdres/VL)*(-PL*(CL - CIlast) + PG*(CG - CL));
        if(CI<0)
        {
        temp = 1/PG+1/PL+1/PPI;
        CI = CG/(PPI*temp);
        CL = CI*(PPI+PL)/PL;
        };
        exponOut[k] = CI*PPI;
        k=k+1;
    }

    for (k=0; k<delaypoint; k++)
        powerLawIn[k] = exponOut[0];
    for (k=delaypoint; k<totalstim*nrep+delaypoint; k++)
        powerLawIn[k] = exponOut[k-delaypoint];
    for (k=totalstim*nrep+delaypoint; k<totalstim*nrep+3*delaypoint; k++)
        powerLawIn[k] = powerLawIn[k-1];

    /* Do lazy downsampling by just taking every tenth sample */
    int idx_down = 0;
    double* sampIHC = (double*)calloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq), sizeof(double));
    for (int i=0; i < (totalstim*nrep+2*delaypoint); i=i+10) {
        idx_down = (int) i/10;
        sampIHC[idx_down] = powerLawIn[i];
    }
    free(powerLawIn); free(exponOut);

    /*----------------------------------------------------------*/
    /*----- Running Power-law Adaptation -----------------------*/
    /*----------------------------------------------------------*/
    k = 0;

    for (indx = 0; indx < floor((totalstim * nrep + 2 * delaypoint) *
                                tdres * sampFreq);
        indx++)
    {
    sout1[k]  = __max( 0, sampIHC[indx] + randNums[indx]- alpha1*I1);
    sout2[k] = __max(0, sampIHC[indx] - alpha2 * I2);

    if (implnt == 1) /* ACTUAL Implementation */
    {
        I1 = 0; I2 = 0;
        for (j=0; j<k+1; ++j)
        {
        I1 += (sout1[j])*binwidth/((k-j)*binwidth + beta1);
        I2 += (sout2[j])*binwidth/((k-j)*binwidth + beta2);
        }
    } /* end of actual */

    if (implnt==0)    /* APPROXIMATE Implementation */
    {
        if (k==0)
        {
        n1[k] = 1.0e-3*sout2[k];
        n2[k] = n1[k]; n3[0]= n2[k];
        }
        else if (k==1)
        {
        n1[k] = 1.992127932802320*n1[k-1]+ 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1]);
        n2[k] = 1.999195329360981*n2[k-1]+ n1[k] - 1.997855276593802*n1[k-1];
        n3[k] = -0.798261718183851*n3[k-1]+ n2[k] + 0.798261718184977*n2[k-1];
        }
        else
        {
        n1[k] = 1.992127932802320*n1[k-1] - 0.992140616993846*n1[k-2]+ 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1] + 0.000000000002347*sout2[k-2]);
        n2[k] = 1.999195329360981*n2[k-1] - 0.999195402928777*n2[k-2]+n1[k] - 1.997855276593802*n1[k-1] + 0.997855827934345*n1[k-2];
        n3[k] =-0.798261718183851*n3[k-1] - 0.199131619873480*n3[k-2]+n2[k] + 0.798261718184977*n2[k-1] + 0.199131619874064*n2[k-2];
        }
        I2 = n3[k];

        if (k==0)
        {
        m1[k] = 0.2*sout1[k];
        m2[k] = m1[k];	m3[k] = m2[k];
        m4[k] = m3[k];	m5[k] = m4[k];
        }
        else if (k==1)
        {
        m1[k] = 0.491115852967412*m1[k-1] + 0.2*(sout1[k] - 0.173492003319319*sout1[k-1]);
        m2[k] = 1.084520302502860*m2[k-1] + m1[k] - 0.803462163297112*m1[k-1];
        m3[k] = 1.588427084535629*m3[k-1] + m2[k] - 1.416084732997016*m2[k-1];
        m4[k] = 1.886287488516458*m4[k-1] + m3[k] - 1.830362725074550*m3[k-1];
        m5[k] = 1.989549282714008*m5[k-1] + m4[k] - 1.983165053215032*m4[k-1];
        }
        else
        {
        m1[k] = 0.491115852967412*m1[k-1] - 0.055050209956838*m1[k-2]+ 0.2*(sout1[k]- 0.173492003319319*sout1[k-1]+ 0.000000172983796*sout1[k-2]);
        m2[k] = 1.084520302502860*m2[k-1] - 0.288760329320566*m2[k-2] + m1[k] - 0.803462163297112*m1[k-1] + 0.154962026341513*m1[k-2];
        m3[k] = 1.588427084535629*m3[k-1] - 0.628138993662508*m3[k-2] + m2[k] - 1.416084732997016*m2[k-1] + 0.496615555008723*m2[k-2];
        m4[k] = 1.886287488516458*m4[k-1] - 0.888972875389923*m4[k-2] + m3[k] - 1.830362725074550*m3[k-1] + 0.836399964176882*m3[k-2];
        m5[k] = 1.989549282714008*m5[k-1] - 0.989558985673023*m5[k-2] + m4[k] - 1.983165053215032*m4[k-1] + 0.983193027347456*m4[k-2];
        }
        I1 = m5[k];
    } /* end of approximate implementation */
    synSampOut[k] = sout1[k] + sout2[k];
    k = k+1;
    }   /* end of all samples */
    free(sout1); free(sout2);
    free(m1); free(m2); free(m3); free(m4); free(m5); free(n1); free(n2); free(n3);

    /*----------------------------------------------------------*/
    /*----- Upsampling to original (High 100 kHz) sampling rate --------*/
    /*----------------------------------------------------------*/
    for(z=0; z<k-1; ++z)
    {
        incr = (synSampOut[z+1]-synSampOut[z])/resamp;
        for(b=0; b<resamp; ++b)
        {
            TmpSyn[z*resamp+b] = synSampOut[z]+ b*incr;
        }
    }
    for (i=0;i<totalstim*nrep;++i)
        synouttmp[i] = TmpSyn[i+delaypoint];

    free(synSampOut); free(TmpSyn); free(sampIHC);
    return((long) ceil(totalstim*nrep));
}