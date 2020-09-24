#ifndef _funzione_
#define _funzione_

#include<cmath>


class FunzioneBase {
	public:
		virtual double Eval(double x) const = 0;
		virtual double EvalSampling(double x) const = 0;
};

class function: public FunzioneBase {

	public:
		double Eval(double x) const { return M_PI*0.5*cos(M_PI*x*0.5); }
		double EvalSampling(double x) const { return ( M_PI*cos(0.5*M_PI*x) ) / (3.* sqrt(1.-x)); }


};

#endif
