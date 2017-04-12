// Platelet-bound chemical species : Purely reaction
#include <algorithm>
#include "ODEsolver.h"
#include "Domain.h"

using namespace std;

namespace Reactions{
	double ODESolver::RK2solve(Domain &domain)
	{

		// total points in domain
		double N = domain.nx * domain.ny;

		// RK-2, need to store k1, k2, k1*dt and k2*dt (temp vars per stage)
		Mesh k1, k2, k1dt, k1plusk2, temp;

		//get the mesh to work on
		Mesh &m = domain.get_mesh();

		// Resize k1 and k2 (HEIGHT x WIDTH)
		k1.resize(domain.ny);
		k2.resize(domain.ny);
		k1dt.resize(domain.ny);
		k1plusk2.resize(domain.ny);
		temp.resize(domain.ny);

		for (int i = 0; i < domain.ny; ++i)
		{
			k1[i].resize(domain.nx);
			k2[i].resize(domain.nx);
			k1dt[i].resize(domain.nx);
			k1plusk2[i].resize(domain.nx);
			temp[i].resize(domain.nx);
		}   

	    double current_time = 0;

	  	// loop over all variables
		while(current_time<T)
		{
			current_time = current_time + dt;
			cout<<"time = "<<current_time<<endl;

			// call function to update bound plateltes z and e put together
			domain.updateTotal();

			//for(int i = 0; i< numV; i++)
			//{
				// for each point in the domain :
				for(int j = 0; j<domain.ny; j++)
					for(int k = 0; k<domain.nx; k++)
					{
						// find k1
						computeK(m[j][k], k1[j][k]);

						// update temp[j][k].Cvals[$] with  m[j][k].Cvals[$] + dt*k1[j][k].Cvals[$]

						// (a) first multiply k1 with dt : vector-scalar multiplication
						std::transform(k1[j][k].Cvals.begin(), k1[j][k].Cvals.end(), k1dt[j][k].Cvals.begin(),
               							std::bind1st(std::multiplies<double>(),dt));

						// (b) http://stackoverflow.com/questions/13728430/element-wise-multiplication-of-two-vectors-in-c
						std::transform( m[j][k].Cvals.begin(), m[j][k].Cvals.end(),
						                k1dt[j][k].Cvals.begin(), std::back_inserter(temp[j][k]. Cvals),
						                std::plus<double>() );
						// find k2
						computeK(temp[j][k], k2[j][k]);

						// apply RK2 formula that uses k1 and k2:
						// (a) k1plusk2 = k1 + k2
						std::transform( k1[j][k].Cvals.begin(), k1[j][k].Cvals.end(),
						                k2[j][k].Cvals.begin(), std::back_inserter(k1plusk2[j][k].Cvals),  
						                std::plus<double>() );
						// (b) multiply (dt/2)*(k1+k2), store in temp
						std::transform(k1plusk2[j][k].Cvals.begin(), k1plusk2[j][k].Cvals.end(), temp[j][k].Cvals.begin(),
               							std::bind1st(std::multiplies<double>(),dt/2));
						// (c) update m[j][k].Cvals[$] with  m[j][k].Cvals[$] + temp
						std::transform( m[j][k].Cvals.begin(), m[j][k].Cvals.end(),
						                temp[j][k].Cvals.begin(), m[j][k].Cvals.begin(),  
						                std::plus<double>() );

					}
			//}

		}
		return current_time;
  
	}

	void ODESolver::computeK(Point &p, Point &kp)
	{

		/*A.9*/ kp.Cvals[PBz2ba] = tchar*(
									consts[k2on]*p.Cvals[PBz2ba]* //k2on*z2*
									(consts[N2b]*Pba + consts[N2se]*Psea - p.Cvals[PBz2ba] - p.Cvals[z2mtot] - p.Cvals[e2mtot])
									- consts[k2off]*p.Cvals[PBz2ba] 
									- consts[kz2promplus]*p.Cvals[PBz2ba]*p.Cvals[PBpro] 
									+ consts[kz2promminus]*p.Cvals[PBz2ba_pro]
								);	
		/*A.10*/kp.Cvals[PBe2ba] = tchar*(
									consts[k2on]*p.Cvals[PBe2ba]* //k2on*e2*
									(consts[N2b]*Pba + consts[N2se]*Psea - p.Cvals[PBz2ba] - p.Cvals[z2mtot] - p.Cvals[e2mtot])
									- consts[k2off]*p.Cvals[PBe2ba] 
									+ consts[kz2promcat]*p.Cvals[PBz2ba_pro]
									+ (consts[kz5e2mcat] + consts[kz5e2mminus])*p.Cvals[PBz5ba_e2ba]
									- consts[kz5e2mplus]*p.Cvals[PBz5ba]*p.Cvals[PBz2ba]
									+ (consts[kz8e2mcat] + consts[kz8e2mminus])*p.Cvals[PBz8ba_e2ba]
									- consts[kz8e2mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe2ba]

								);
		/*A.11*/kp.Cvals[PBz5ba] = tchar*(
									consts[k5on]*p.Cvals[PBz5ba]* //k5on*z5*
									(consts[N5b]*Pba + consts[N5se]*Psea - p.Cvals[ze5mtot])
									- consts[k5off]*p.Cvals[PBz5ba]
									- consts[kz5e10mplus]*p.Cvals[PBz5ba]*p.Cvals[PBe10ba]
									+ consts[kz5e10mminus]*p.Cvals[PBz5ba_e10ba]
									- consts[kz5e2mplus]*p.Cvals[PBz5ba]*p.Cvals[PBe2ba]
									+ consts[kz5e2mminus]*p.Cvals[PBz5ba_e2ba]

								);
		/*A.12*/kp.Cvals[PBe5ba] = tchar*(
									consts[k5on]*p.Cvals[PBe5ba]* //k5on*e5*
									(consts[N5b]*Pba + consts[N5se]*Psea - p.Cvals[ze5mtot])
									- consts[k5off]*p.Cvals[PBe5ba]
									+ consts[kz5e10mcat]*p.Cvals[PBz5ba_e10ba]
									+ consts[kz5e2mcat]*p.Cvals[PBz5ba_e2ba]
									+ consts[kprominus]*p.Cvals[PBpro]
									- consts[kproplus]*p.Cvals[PBe5ba]*p.Cvals[PBe10ba]
									- consts[kapce5mplus]*p.Cvals[PBapc]*p.Cvals[PBe5ba]
									+ consts[kapce5mminus]*p.Cvals[PBapc_e5ba]
								);
		/*A.13*/kp.Cvals[PBz8ba] = tchar*(
									consts[k8on]*p.Cvals[PBz8ba]* //k8on*z8*
									(consts[N8b]*Pba + consts[N8se]*Psea - p.Cvals[ze8mtot])
									- consts[k8off]*p.Cvals[PBz8ba]
									- consts[kz8e10mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe10ba]
									+ consts[kz8e10mminus]*p.Cvals[PBz8ba_e10ba]
									- consts[kz8e2mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe2ba]
									+ consts[kz8e2mminus]*p.Cvals[PBz8ba_e2ba]
								);
		/*A.14*/kp.Cvals[PBe8ba] = tchar*(
									consts[k8on]*p.Cvals[PBe8ba]* //k8on*e8*
									(consts[N8b]*Pba + consts[N8se]*Psea - p.Cvals[ze8mtot])
									- consts[k8off]*p.Cvals[PBe8ba]
									+ consts[kz8e10mcat]*p.Cvals[PBz8ba_e10ba]
									+ consts[kz8e2mcat]*p.Cvals[PBz8ba_e2ba]
									+ consts[ktenminus]*p.Cvals[PBten]
									- consts[ktenplus]*p.Cvals[PBe8ba]*p.Cvals[PBe9ba]
									- consts[kapce8mplus]*p.Cvals[PBapc]*p.Cvals[PBe8ba]
									+ consts[kapce8mminus]*p.Cvals[PBapc_e8ba]
								);
		/*A.15*/kp.Cvals[PBz9ba] = tchar*(
									consts[k9on]*p.Cvals[PBz9ba]* //k9on*z9*
									(consts[N9b]*Pba + consts[N9se]*Psea - p.Cvals[ze9mtot])
									-consts[k9off]*p.Cvals[PBz9ba]
								);
		/*A.16*/kp.Cvals[PBe9ba] = tchar*(
									consts[k9on]*p.Cvals[PBe9ba]* // k9on*e9*
									(consts[N9b]*Pba + consts[N9se]*Psea - p.Cvals[ze9mtot])
									-consts[k9off]*p.Cvals[PBe9ba]
									+ consts[ktenminus]*p.Cvals[PBten] 
									- consts[ktenplus]*p.Cvals[PBe8ba]*p.Cvals[PBe9ba]
								);
		/*A.17*/kp.Cvals[PBz10ba] = tchar*(
									consts[k10on]*p.Cvals[PBz10ba]* //k10on*z10*
									(consts[N10b]*Pba + consts[N10se]*Psea - p.Cvals[ze10mtot])
									- consts[k10off]*p.Cvals[PBz10ba]
									+ consts[kz10tenmminus]*p.Cvals[PBz10ba_ten]
									- consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBten]
									+ consts[kz10tenmminus]*p.Cvals[PBz10ba_tenstar]
									- consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBtenstar]
								);
		/*A.18*/kp.Cvals[PBe10ba] = tchar*(
									consts[k10on]*p.Cvals[PBe10ba]* // k10on*e10*
									(consts[N10b]*Pba + consts[N10se]*Psea - p.Cvals[ze10mtot])
									- consts[k10off]*p.Cvals[PBe10ba]
									+ (consts[kz5e10mminus] + consts[kz5e10mcat])*p.Cvals[PBz5ba_e10ba]
									- consts[kz5e10mminus]*p.Cvals[PBz5ba]*p.Cvals[PBe10ba]
									+ (consts[kz8e10mminus] + consts[kz8e10mcat])*p.Cvals[PBz8ba_e10ba]
									- consts[kz8e10mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe10ba]
									+ consts[kprominus]*p.Cvals[PBpro]
									- consts[kproplus]*p.Cvals[PBe5ba]*p.Cvals[PBe10ba]
									+ consts[kz10tenmcat]*p.Cvals[PBz10ba_ten]
									+ consts[kz10tenmcat]*p.Cvals[PBz10ba_tenstar]
								);
		/*A.19*/kp.Cvals[PBten] = tchar*(
									consts[ktenplus]*p.Cvals[PBe8ba]*p.Cvals[PBe9ba]
									- consts[ktenminus]*p.Cvals[PBten]
									- consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBten]
									+ (consts[kz10tenmcat] + consts[kz10tenmminus])*p.Cvals[PBz10ba_ten]
								);
		/*A.20*/kp.Cvals[PBpro] = tchar*(
									consts[kproplus]*p.Cvals[PBe5ba]*p.Cvals[PBe10ba]
									- consts[kprominus]*p.Cvals[PBpro]
									- consts[kz2promplus]*p.Cvals[PBz2ba]*p.Cvals[PBpro]
									+ (consts[kz2promcat] + consts[kz2promminus])*p.Cvals[PBz2ba_pro]
								);
		/*A.21*/kp.Cvals[PBz2ba_pro] = tchar*(
									consts[kz2promplus]*p.Cvals[PBz2ba]*p.Cvals[PBpro]
									- (consts[kz2promminus] + consts[kz2promcat])*p.Cvals[PBz2ba_pro]
								);
		/*A.22*/kp.Cvals[PBz5ba_e2ba] = tchar*(
									consts[kz5e2mplus]*p.Cvals[PBz5ba]*p.Cvals[PBe2ba]
									- (consts[kz5e2mminus] + consts[kz5e2mcat])*p.Cvals[PBz5ba_e2ba]
								);
		/*A.23*/kp.Cvals[PBz5ba_e10ba] = tchar*(
									consts[kz5e10mplus]*p.Cvals[PBz5ba]*p.Cvals[PBe10ba]
									- (consts[kz5e10mminus] + consts[kz5e10mcat])*p.Cvals[PBz5ba_e10ba]
								);
		/*A.24*/kp.Cvals[PBz8ba_e2ba] = tchar*(
									consts[kz8e2mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe2ba]
									- (consts[kz8e2mminus] + consts[kz8e2mcat])*p.Cvals[PBz8ba_e2ba]
								);
		/*A.25*/kp.Cvals[PBz8ba_e10ba] = tchar*(
									consts[kz8e10mplus]*p.Cvals[PBz8ba]*p.Cvals[PBe10ba]
									- (consts[kz8e10mminus] + consts[kz8e10mcat])*p.Cvals[PBz8ba_e10ba]
								);
		/*A.26*/kp.Cvals[PBz10ba_ten] = tchar*(
									consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBten]
									- (consts[kz10tenmminus] + consts[kz10tenmcat])*p.Cvals[PBz10ba_ten]
								);
		/*A.27*/kp.Cvals[PBapc_e5ba] = tchar*(
									consts[kapce5mplus]*p.Cvals[PBapc]*p.Cvals[PBe5ba] - (consts[kapce5mminus] + consts[kapce5mcat])*p.Cvals[PBapc_e5ba]
								);
	//	/*A.27*/kp.Cvals[PBapc_e5ba] = tchar*(
  	//								Consts::kapce5mplus*p.PBapc*p.PBe5ba - (Consts::kapce5mminus + Consts::kapce5mcat)*p.PBapc_e5ba
	//							);
		/*A.28*/kp.Cvals[PBapc_e8ba] = tchar*(
									consts[kapce8mplus]*p.Cvals[PBapc]*p.Cvals[PBe8ba]
									- (consts[kapce8mminus] + consts[kapce8mcat])*p.Cvals[PBapc_e8ba]
								);
		/*A.29*/kp.Cvals[PBe9starba] = tchar*(
									consts[k9on]*p.Cvals[PBe9ba]* //k9on*e9*
									(consts[N9starb]*Pba + consts[N9starse]*Psea - p.Cvals[PBe9starba] - p.Cvals[PBtenstar] - p.Cvals[PBz10ba_tenstar])
									- consts[k9off]*p.Cvals[PBe9starba]
									+ consts[ktenminus]*p.Cvals[PBtenstar]
									- consts[ktenplus]*p.Cvals[PBe8ba]*p.Cvals[PBe9ba]
								);
		/*A.30*/kp.Cvals[PBtenstar] = tchar*(
									consts[ktenplus]*p.Cvals[PBe8ba]*p.Cvals[PBe9starba]
									- consts[ktenminus]*p.Cvals[PBtenstar]
									+ (consts[kz10tenmminus] + consts[kz10tenmcat])*p.Cvals[PBz10ba_tenstar]
									- consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBtenstar]
								);
		/*A.31*/kp.Cvals[PBz10ba_tenstar] = tchar*(
									consts[kz10tenmplus]*p.Cvals[PBz10ba]*p.Cvals[PBtenstar] 
									- (consts[kz10tenmminus] + consts[kz10tenmcat])*p.Cvals[PBz10ba_tenstar]
								);
	}
}


