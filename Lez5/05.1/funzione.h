#ifndef _funzione_
#define _funzione_

#include<cmath>

double a0 = 0.0529E-9; //Bohr Radius

class FunzioneBase {
	public:
		virtual double WF_1s(double x, double y, double z) const = 0;
		virtual double WF_2p(double x, double y, double z) const = 0;
};

class WF: public FunzioneBase {  //Using Bohr Radius Units

	public:   //Normalization factors are not considered
		double WF_1s (double x, double y, double z) const { return exp( -sqrt(x*x+y*y+z*z) ); }
		double WF_2p (double x, double y, double z) const { return sqrt(x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z)/2.)*cos(atan(sqrt(x*x+y*y)/z)); }

};

#endif
