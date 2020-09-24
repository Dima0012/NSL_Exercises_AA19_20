#include"integraleMC.h"

double IntegraleMC::IntegraleHoM(double xmin, double xmax, double fmax, FunzioneBase *f, int punti){

	int Nhit = 0;	

	for(int i = 0; i<punti; i++){

		double x = m_rnd.Rannyu(xmin, xmax);
		double y = m_rnd.Rannyu(0, fmax);

		if (y < f->Eval(x)) Nhit++;

	}

	return (xmax-xmin)*fmax*(double(Nhit)/double(punti));

}

double IntegraleMC::IntegraleAVE(double xmin, double xmax, FunzioneBase *f, int punti){

	double ftot = 0.;

	for(int i = 0; i<punti; i++){
		ftot += f->Eval(m_rnd.Rannyu(xmin, xmax));
	}
	
	if(punti == 0) return 0;

	return (xmax-xmin)*ftot/punti;

}

double IntegraleMC::IntegraleSAMP(double xmin, double xmax, FunzioneBase *f, int punti){

	double ftot = 0.;

	for(int i = 0; i<punti; i++){

		//double x = 1. - pow( m_rnd.Rannyu(0., 1.) + 1. , 2./3.); 	//Distribuzione sampling (???)

		double x = 0.;
		double r = 0.; 

		while (true) {
			x = m_rnd.Rannyu(xmin, xmax);
			r = m_rnd.Rannyu();

			double val = 1.5*sqrt(1.-x);
			double pmax = 1.5;

			if ( r < val/pmax )  	//Accept-Reject
				break;
		}

		ftot += f->EvalSampling(x);
	}

	return (xmax-xmin)*ftot/punti;
}


