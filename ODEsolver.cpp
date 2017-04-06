// Platelet-bound chemical specied : Purely reaction

#include "ODEsolver.h"
#include "Domain.h"

using namespace std;

namespace ODESolver{

// signature of function to display output concentrations:
void displayC(ODESolver *, Domain*, double );

/////////////////////////////
int main(int argc, char *argv[])
{
	ODESolver *solver = new ODESolver(2,0.5);

	Domain *domain = new Domain(50,20,100,80);

	// total points in domain
	double N = domain->nx * domain->ny;
	double Pba, Psea;
	double current_time = 0;

	ODESolver::ConstMap &consts = solver->getConsts();

	// check assignment of constants
	/*cout<<"Constants for ODEs assigned as follows:"<<endl;
	for (const auto &kv : consts)
	{
			cout << kv.first<<" = "<< kv.second <<" " <<consts[ODESolver::k2on]<< endl;
	}*/

	cout<< "Enter Pba and Psea:"<<endl;
	cin>>Pba>>Psea;
	cout<< "You entered: Pba="<<Pba<< " and Psea="<<Psea<<endl;;

	// RK-2, need to store k1 and k2 (temp vars per stage)
	std::vector<std::vector<Point>> k1;
	std::vector<std::vector<Point>> k2;

	// Resize k1 and k2 (HEIGHT x WIDTH)
	k1.resize(domain->ny);
	k2.resize(domain->ny);
	for (int i = 0; i < domain->ny; ++i)
	{
		k1[i].resize(domain->nx);
		k2[i].resize(domain->nx);
	}   

	// begin time tracking:
	time_t start = time(NULL);

	// loop over all variables
	while(current_time<solver->T)
	{
		current_time = current_time + solver->dt;

		// call function to update bound plateltes z and e put together
		updateTotal();

		for( int i = 0; i< solver->numV; i++)
		{
			// for each point in the domain :
			for(int j = 0; j<domain->ny; j++)
				for(int k = 0; k<domain->nx;k++)
						k1[j][k].Cvals[ODESolver::PBz2ba] = solver->dt * 
						(
								consts[ODESolver::k2on]*domain->mesh[j][k].Cvals[ODESolver::PBz2ba]* //k21on*z2*
								(consts[ODESolver::N2b]*Pba + consts[ODESolver::N2se]*Psea - domain->mesh[j][k].Cvals[ODESolver::PBz2ba] - domain->mesh[j][k].Cvals[ODESolver::z2mtot] - domain->mesh[j][k].Cvals[ODESolver::e2mtot])
						);
		}

	}

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

}