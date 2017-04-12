// Platelet-bound chemical species : Purely reaction

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
			mesh[i].resize(domain.nx);
		}   

	    double current_time = 0;

	  	// loop over all variables
		while(current_time<T)
		{
			current_time = current_time + dt;

			// call function to update bound plateltes z and e put together
			domain->updateTotal();

			for(int i = 0; i< numV; i++)
			{
				// for each point in the domain :
				for(int j = 0; j<domain.ny; j++)
					for(int k = 0; k<domain.nx; k++)
					{
						// find k1
						computeK(m[j][k], k1[j][k]);

						// update temp[j][k].Cvals[$] with  m[j][k].Cvals[$] + dt*k1[j][k].Cvals[$]

						// (a) first multiply k1 with dt : vector-scalar multiplication
						std::transform(k1[j][k].Cvals.begin(), k1[j][k].Cvals.end(), k1dt[j][k].Cvals.begin(),
               							std::bind1st(std::multiplies<T>(),dt));

						// (b) http://stackoverflow.com/questions/13728430/element-wise-multiplication-of-two-vectors-in-c
						std::transform( m[j][k].Cvals.begin(), m[j][k].Cvals.end(),
						                k1dt[j][k].begin(), temp[j][k].Cvals.begin(),  
						                std::adds<int>() );
						// find k2
						computeK(temp[j][k], k2[j][k]);

						// apply RK2 formula that uses k1 and k2:
						// (a) k1plusk2 = k1 + k2
						std::transform( k1[j][k].Cvals.begin(), k1[j][k].Cvals.end(),
						                k2[j][k].begin(), k1plusk2[j][k].Cvals.begin(),  
						                std::adds<int>() );
						// (b) multiply (dt/2)*(k1+k2), store in temp
						std::transform(k1plusk2[j][k].Cvals.begin(), k1plusk2[j][k].Cvals.end(), temp[j][k].Cvals.begin(),
               							std::bind1st(std::multiplies<T>(),dt/2));
						// (c) update m[j][k].Cvals[$] with  m[j][k].Cvals[$] + temp
						std::transform( m[j][k].Cvals.begin(), m[j][k].Cvals.end(),
						                temp[j][k].begin(), m[j][k].Cvals.begin(),  
						                std::adds<int>() );

					}
			}

		}
		return current_time;
  
	}

	void computeK(Point &p, Point &kp)
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
									
								);
	}
}

//move this to another file: driver.cpp

using namespace Reactions;

// signature of function to display output concentrations:
void displayC(ODESolver *, Domain*, double );

/////////////////////////////
int main(int argc, char *argv[])
{
	ODESolver *solver = new ODESolver(2,0.5);

	Domain domain = new Domain(50,20,100,80);
	double current_time = 0; 

	double Pba, Psea;

	ConstMap &consts = solver->getConsts();

	// check assignment of constants
	/*cout<<"Constants for ODEs assigned as follows:"<<endl;
	for (const auto &kv : consts)
	{
			cout << kv.first<<" = "<< kv.second <<" " <<consts[ODESolver::k2on]<< endl;
	}*/

	cout<< "Enter Pba and Psea:"<<endl;
	cin>>Pba>>Psea;
	cout<< "You entered: Pba="<<Pba<< " and Psea="<<Psea<<endl;;

	// set Pba and Psea for solver
	solver->setPba(Pba);
	solver->setPsea(Psea);

	// begin time tracking:
	time_t start = time(NULL);

	current_time = solver->RK2solve(domain);

	cout<<"time elapsed="<<time(NULL) - start<<endl;;

	displayC(solver, domain, current_time);


	return 0;
}

void displayC(ODESolver  *solver, Domain *domain, double current_time){
	cout<<"\nComputation complete. Here are the concentrations for chemicals at a point (x,y) in time "
	<<current_time<<" s"<<endl;
	cout<<" PBz2ba, PBe2ba, PBz5ba, PBe5ba, PBz8ba, PBe8ba, PBz9ba, PBe9ba,PBz10ba, PBe10ba, PBten,\
	PBpro, PBz2ba_pro, PBz5ba_e2ba, PBz5ba_e10ba, PBz8ba_e2ba, PBz8ba_e10ba, PBz10ba_ten, PBapc_e5ba,\
	PBapc_e8ba, PBe9starba, PBtenstar, PBz10ba_tenstar, PBPltsADP, PBPltse2, e2mtot, z2mtot, z5mtot, \
	e5mtot, z8mtot, e8mtot, z9mtot, e9mtot, e10mtot, z10mtot :"<<endl;
	for( int i = 0; i< solver->numV; i++)
	{
		cout<<domain->mesh[2][2].Cvals[i]<<", ";
	}

}

