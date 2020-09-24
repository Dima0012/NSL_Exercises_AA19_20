#ifndef _integraleMC_
#define _integraleMC_

#include"random.h"
#include"funzione.h"

class IntegraleMC {

	public:
		IntegraleMC(Random rnd) {m_rnd = rnd;} 

		double IntegraleHoM(double xmin, double xmax, double fmax, FunzioneBase *f, int punti);
		double IntegraleAVE(double xmin, double xmax, FunzioneBase *f, int punti);
		double IntegraleSAMP(double xmin, double xmax, FunzioneBase *f, int punti);
	
	private:
		Random m_rnd;

};

#endif
