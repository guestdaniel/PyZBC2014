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


#include <stdlib.h>
#include <math.h>
#include "complex.hpp"

/* divide */
COMPLEX compdiv(COMPLEX ne,COMPLEX de)
{
  double d;
  COMPLEX z;
  d=de.x*de.x+de.y*de.y;
  z.x=(ne.x*de.x+ne.y*de.y)/d;
  z.y=(ne.y*de.x-ne.x*de.y)/d;
  return(z);
}
/* this returns a complex number equal to exp(i*theta) */
COMPLEX compexp(double theta)
{
  COMPLEX dummy;
  dummy.x = cos(theta);
  dummy.y = sin(theta);
  return dummy;
}
/* Multiply a complex number by a scalar */
COMPLEX compmult(double scalar, COMPLEX compnum)
{
 COMPLEX answer;
 answer.x = scalar * compnum.x;
 answer.y = scalar * compnum.y;
 return answer;
}
/* Find the product of 2 complex numbers */
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2)
{
 COMPLEX answer;
 answer.x = (compnum1.x * compnum2.x) - (compnum1.y * compnum2.y);
 answer.y = (compnum1.x * compnum2.y) + (compnum1.y * compnum2.x);
 return answer;
}
/* add 2 complex numbers */
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x;
 answer.y = summand1.y + summand2.y;
 return answer;
}
/* add three complex numbers */
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x + summand3.x;
 answer.y = summand1.y + summand2.y + summand3.y;
 return answer;
}

/* subtraction: complexA - complexB */
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB)
{
 COMPLEX answer;
 answer.x = complexA.x - complexB.x;
 answer.y = complexA.y - complexB.y;
 return answer;
}
/* Get the real part of the complex */
double REAL(COMPLEX compnum)
{ return(compnum.x); } 

/* Get the imaginary part of the complex */
double IMAG(COMPLEX compnum)
{ return(compnum.y); } 

/* Get the conjugate of the complex signal */
COMPLEX compconj(COMPLEX complexA)
{
  COMPLEX answer;
  answer.x = complexA.x;
  answer.y = -complexA.y;
  return (answer);
}
