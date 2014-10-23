/*   Class interface definition for func_kolmogorov.cpp
 *   VERSION 0.1
 *
 *   A class derived from FunctionObject (function_object.h),
 * which produces the luminosity as a function of radius for an elliptical
 * Kolmogorov profile, approximated by the sum of two Moffat functions.
 *
 * PARAMETERS:
 * PA = params[0 + offsetIndex];
 * ell = params[1 + offsetIndex];
 * I_0 = params[2 + offsetIndex ];
 * fwhm = params[3 + offsetIndex ];
 *
 *
 */


// CLASS Kolmogorov:

#include "function_object.h"

//#define CLASS_SHORT_NAME  "Kolmogorov"


class Kolmogorov : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    Kolmogorov( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r );
    int  CalculateSubsamples( double r );


  private:
    double  x0, y0, PA, ell, I_0, fwhm;   // parameters
    double  alpha7, alpha2;
    double  q, PA_rad, cosPA, sinPA;   // other useful (shape-related) quantities
};
