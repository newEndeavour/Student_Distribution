/*
  File:         Student_Distribution.cpp
  Version:      0.0.2
  Date:         23-Jan-2019
  Revision:     30-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to Student_Distribution.h for Edition Notes.

  Student_Distribution.cpp - Library for 'duino
  https://github.com/newEndeavour/Student_Distribution

  Student_Distribution implements a Student distribution. 

  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "Arduino.h"
#include <math.h>
#include <float.h>     // required for DBL_MAX and LDBL_MAX
#include "Student_Distribution.h"
#include "Hypergeometric_Function.h"
#include "Beta_Function.h"
#include "Gamma_Function.h"


// Constructor /////////////////////////////////////////////////////////////////
// Function that handles the creation and setup of instances

Student_Distribution::Student_Distribution(double _Nu)
{

	// Object parameter's error handling
	error = 1;

	if (_Nu<=0) 	error = -2;

	//Set initial values	
	Nu			= _Nu;			// 

}


// Public Methods //////////////////////////////////////////////////////////////
//Probability Density Function
double Student_Distribution::GetPDF(double x)
{
double ln_density;

	if (error<0)
		return error;

	//Particular cases
	if (Nu==1) 
		return 1/(CONSTANT_Pi * (1+x*x));

	if (Nu==2) 
		return 1/(pow(1+x*x,3/2));

	//General cases
	//return Gamma_Function(.5*(Nu+1)) * pow(1+pow(x,2)/Nu,-(Nu+1)*0.5) / (sqrt(Nu*CONSTANT_Pi) * Gamma_Function(.5*Nu));

	ln_density = -(double)(Nu+1)/2.0 * log(1.0 + x * x /(double)Nu)
        	     - 0.5*log((double)Nu)
                     - Ln_Beta_Function(0.5 * (double)Nu, 0.5);

	return exp(ln_density);
}


//Cumulative Distribution Function
double Student_Distribution::GetCDF(double x)
{
double a = (double) Nu / 2.0;
double beta;

	if (error<0)
		return error;

	//Particular cases
	if (Nu==1)
		return 0.5 + atan(x) / CONSTANT_Pi;			//atan() is the arc tangent function
	if (Nu==2) 
		return 0.5 + x / (2.0 * sqrt(2+x*x));

	//General cases
	//F(x) = 1/2 + x.G((Nu+1)/2 * 1F2(1/2,(Nu+1)/2,3/2,-x^2/Nu) / (sqrt(Nu*Pi) * G(Nu/2))
	//return ( 0.5 + x*Gamma_Function(.5*(Nu+1)) * 
	//	HyperGeometric(.5,.5*(Nu+1),1.5,-x*x/Nu)/(sqrt(Nu*CONSTANT_Pi) * Gamma_Function(.5*Nu)) );
	
	beta = Beta_Distribution(1/(1+x*x/Nu),a,0.5);
	if (x>0.0) 
		return 1-.5*beta;
   	else if (x<0.0) 
		return .5*beta;
   	
	return 0.5;
}


//Mean
double 	Student_Distribution::GetMean(void)
{
	if (error<0)
		return error;

	//0 for ν > 1, otherwise undefined
	if (Nu > 1) 
		return 0.0;
	else
		return -9999;

}


//Variance
double 	Student_Distribution::GetVariance(void)
{
	if (error<0)
		return error;

	// ν / (ν - 2) if nu > 2, otherwise undefined
	if (Nu > 2) 
		return Nu / (Nu - 2);
	else
		if (Nu>1)
			return DBL_MAX;
		else
			return -9999;

}


//Std Deviation
double 	Student_Distribution::GetStdDeviation(void)
{
	if (error<0)
		return error;

	// ν / (ν - 2) if nu > 2, otherwise undefined
	if (Nu > 2) 
		return sqrt(Nu / (Nu - 2));
	else 
		if (Nu>1)
			return DBL_MAX;
		else
			return -9999;
}


//Skewness
double 	Student_Distribution::GetSkewness(void)
{
	if (error<0)
		return error;

	//0 for ν > 3, otherwise undefined
	if (Nu > 3) 
		return 0.0;
	else
		return -9999;
}


//Kurtosis
double 	Student_Distribution::GetKurtosis(void)
{
	if (error<0)
		return error;

	// 6 / (ν - 4) if nu > 4, otherwise undefined
	if (Nu > 4) 
		return 6.0 / (Nu - 4);
	else 
		if (Nu > 2)
			return DBL_MAX;
		else
			return -9999;
}


//Arc Tan function
double Student_Distribution::GetArcTan(double x)
{
	return atan(x);
}


//Return Quantile z(P) from probability P
double Student_Distribution::GetQuantile(double p)
{
double Vm;
double Vh = 300;
double Vl = -300;
double Pr;
int i = 0;
double Eps;
	
	if (p <= 0.0) {
		return Vl;
	} else if (p >= 1.0) {
        	return Vh;
	} else {        
        	do 
		{
          		i++;
          		Vm = (Vh + Vl) / 2;
            
			Pr = GetCDF(Vm);
          		Eps = abs(Pr - p);
        
          		//New Boundary selection
          		if (Pr > p) {
				Vh = Vm;
          		} else {
				Vl = Vm;
			}
            
        	} 
		while ((Eps > CONSTANT_EpsStop) && (i < 70));
	}
            
        if (i >= 70) {
            return -9999;
        } else {
            return Vm;
    	}

}


double Student_Distribution::GetNu(void)
{
	return Nu;
}


// Private Methods /////////////////////////////////////////////////////////////
// Functions only available to other functions in this library


// /////////////////////////////////////////////////////////////////////////////

