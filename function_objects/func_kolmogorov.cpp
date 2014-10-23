/* FILE: func_moffat.cpp ----------------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for a Kolmogorov PSF modeled as a sum of
 *   two Moffat functions: one with beta = 7 and I_0 = 0.86, and one
 *   with beta = 2 and I_0 = 0.14.
 *   This will generate a normalized profile (total flux = 1.0).
 *
 *   NOTES: This function returns the intensity using the Moffat function:
 *      I(r) = I_0 / [1 + (r/alpha)^2]^beta
 *
 *   User inputs are I_0 and FWHM; For the specific Moffat notes, see
 *   the comment section in func_moffat.cpp.
 *
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x and y.
 *      GetValue() then completes the calculation, using the actual value
 *      of x and y, and returns the result.
 *      So for an image, we expect the user to call Setup() once at
 *      the start, then loop through the pixels of the image, calling
 *      GetValue() to compute the function results for each pixel coordinate
 *      (x,y).
 *
 *   NOTE: Currently, we assume input PA is in *degrees* [and then we
 * convert it to radians] relative to +x axis.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]: 23 Oct 2014: Created (as modification of func_moffat.cpp).
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>

#include "func_kolmogorov.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 4;
const char  PARAM_LABELS[][20] = {"PA", "ell", "I_0", "fwhm"};
const char  FUNCTION_NAME[] = "Kolmogorov PSF";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char Kolmogorov::className[] = "Kolmogorov";

/* ---------------- Useful macros for Moffat function ------------------ */

double inline get_alpha(double fwhm, double beta)
{
  double exponent = pow(2.0, 1.0/beta);
  return 0.5*fwhm/sqrt(exponent - 1.0);
}

double inline moffat(double r, double I_0, double alpha, double beta)
{
  double  scaledR, denominator;

  scaledR = r / alpha;
  denominator = pow((1.0 + scaledR*scaledR), beta);
  return (I_0 / denominator);
}

/* ---------------- CONSTRUCTOR ---------------------------------------- */

Kolmogorov::Kolmogorov( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void Kolmogorov::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  ell = params[1 + offsetIndex];
  I_0 = params[2 + offsetIndex];
  fwhm = params[3 + offsetIndex];

  // pre-compute useful things for this round of invoking the function
  q = 1.0 - ell;
  // convert PA to +x-axis reference and then to radians
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  // compute alpha:
  alpha7 = get_alpha(fwhm, 7.0);
  alpha2 = get_alpha(fwhm, 2.0);
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for a Kolmogorov PSF function at radius r,
// with the various parameters and derived values (alpha, PA, etc.)
// pre-calculated by Setup().

double Kolmogorov::CalculateIntensity( double r )
{
  double moffat7, moffat2;

  moffat7 = moffat(r, 0.86 * I_0, alpha7, 7.0);
  moffat2 = moffat(r, 0.14 * I_0, alpha2, 2.0);
  return moffat7 + moffat2;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double Kolmogorov::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp_scaled, r, totalIntensity;
  int  nSubsamples;
  
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio
  xp = x_diff*cosPA + y_diff*sinPA;
  yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
  r = sqrt(xp*xp + yp_scaled*yp_scaled);
  
  nSubsamples = CalculateSubsamples(r);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-picel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        xp = x_diff*cosPA + y_diff*sinPA;
        yp_scaled = (-x_diff*sinPA + y_diff*cosPA)/q;
        r = sqrt(xp*xp + yp_scaled*yp_scaled);
        theSum += CalculateIntensity(r);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(r);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// Moffat function.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int Kolmogorov::CalculateSubsamples( double r )
{
  int  nSamples = 1;
  
  // Chien Peng algorithm for Moffat function
  if ((doSubsampling) && (r < 10.0)) {
    if ((alpha7 <= 1.0) && (r <= 1.0))
      nSamples = min(100, (int)(4 * SUBSAMPLE_R / alpha7));
    else {
      if (r <= 3.0)
        nSamples = 2 * SUBSAMPLE_R;
      else
        nSamples = min(100, (int)(2 * SUBSAMPLE_R / r));
    }
  }
  return nSamples;
}



/* END OF FILE: func_kolmogorov.cpp ---------------------------------------- */
