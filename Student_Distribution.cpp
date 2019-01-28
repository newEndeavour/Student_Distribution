/*
  File:         Student_Distribution.cpp
  Version:      0.0.1
  Date:         23-Jan-2019
  Revision:     23-Jan-2019
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
#include "Student_Distribution.h"
#include "Gamma_Function.h"


// Constructor /////////////////////////////////////////////////////////////////
// Function that handles the creation and setup of instances

Student_Distribution::Student_Distribution(float _Nu)
{

	// Object parameter's error handling
	error = 1;

	//Set initial values	
	Nu			= _Nu;			// 

}


// Public Methods //////////////////////////////////////////////////////////////
//Probability Density Function
float Student_Distribution::GetPDF(float x)
{
    return ( Gamma_Function(Nu/2+1/2) * (pow(1+pow(x,2)/Nu,-(Nu+1)/2)) ) / (sqrt(Nu*CONSTANT_Pi) * Gamma_Function(Nu/2));
}


//TEMPORARY
double HyperGeometric_Function(double a, double b, double c, double d) 
{
	return 1;
}


//Cumulative Distribution Function
//
//	F(x) = 1/2 + x.G((Nu+1)/2 * 1F2(1/2,(Nu+1)/2,3/2,-x^2/Nu) / (sqrt(Nu*Pi) * G(Nu/2))
//
// 	Where:
//	- 1F2 	: Hypergeometric function
//	- G	: Gamma function
//	 
float Student_Distribution::GetCDF(float x)
{
    return ( 1/2 + x*Gamma_Function(Nu/2+1/2) * 
	HyperGeometric_Function(1/2,Nu/2+1/2,3/2,-pow(x,2)/Nu)/(sqrt(Nu*CONSTANT_Pi) * Gamma_Function(Nu/2)) );
}


//Return Quantile z(P) from probability P
float Student_Distribution::GetQuantile(float p)
{
float Vm;
float Vh = 16;
float Vl = -16;
float Pr;
int i = 0;
float Eps;

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
			
			/*
			//DEBUG
			Serial.print("\nF(x) | ");
			Serial.print("i:");
			Serial.print(i);
			Serial.print("\tVl:");
			Serial.print(Vl,4);
			Serial.print("\tVh:");
			Serial.print(Vh,4);
			Serial.print("\tVm:");
			Serial.print(Vm,4);
			Serial.print("\t ->F(Vm):");
			Serial.print(Pr,4);
			Serial.print("\t ->Eps:");
			Serial.print(Eps,4);
			*/
        
          		//New Boundary selection
          		if (Pr > p) {
				Vh = Vm;
				//Serial.print("\t (Pr > p) Vh=Vm->Vh:");
				//Serial.print(Vh,4);
          		} else {
				Vl = Vm;
				//Serial.print("\t (Pr < p) Vl=Vm->Vl:");
				//Serial.print(Vl,4);
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


float Student_Distribution::GetNu(void)
{
	return Nu;
}


// Private Methods /////////////////////////////////////////////////////////////
// Functions only available to other functions in this library


// /////////////////////////////////////////////////////////////////////////////

