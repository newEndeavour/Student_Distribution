/*
  File:         Student_Distribution.h
  Version:      0.0.2
  Date:         27-Jan-2019
  Revision:     29-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Student_Distribution.h - Library for 'duino
  https://github.com/newEndeavour/Student_Distribution

  Student_Distribution implements the Student T distribution 
  with v (Nu) degrees of Freedom.

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

  Editions:
  - 0.0.1	: First version
  - 0.0.2	: Corrections for PDF and CDF for particular cases where Nu=1 (the
		  Cauchy distribution) or Nu=2. 
		  Added ArcTan function.
		  Extended the validity domain of boundaries for GetQuantile.
		  

*/


// ensure this library description is only included once
#ifndef Student_Distribution_h
#define Student_Distribution_h


#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif


//Student Constant Parameters
#define CONSTANT_EpsStop 0.0000001
    
//Number Pi
#define CONSTANT_Pi 3.14159265358979


// library interface description
class Student_Distribution
{
  // user-accessible "public" interface
  public:
  // methods
	Student_Distribution(double _Nu);
	
	double 	GetPDF(double x);
	double 	GetCDF(double x);

	double 	GetMean(void);
	double 	GetVariance(void);
	double 	GetStdDeviation(void);
	double 	GetSkewness(void);
	double 	GetKurtosis(void);

	double 	GetQuantile(double p);

	double 	GetArcTan(double x);
	double 	GetNu(void);

  // library-accessible "private" interface
  private:
  // variables
	int 	error;

	double	Nu;			// 
	
  // methods
};

#endif
