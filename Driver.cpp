#include <iostream>
#include "ODEsolver.h"
#include "Domain.h"

using namespace std;
using namespace Reactions;

// signature of function to display output concentrations:
void displayC(ODESolver *, Domain&, double );

/////////////////////////////
int main(int argc, char *argv[])
{
	ODESolver *solver = new ODESolver(2,0.5);

	Domain domain(50,20,100,80);
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

void displayC(ODESolver  *solver, Domain &domain, double current_time){
	cout<<"\nComputation complete. Here are the concentrations for chemicals at a point (x,y) in time "
	<<current_time<<" s"<<endl;
	cout<<" PBz2ba, PBe2ba, PBz5ba, PBe5ba, PBz8ba, PBe8ba, PBz9ba, PBe9ba,PBz10ba, PBe10ba, PBten,\
	PBpro, PBz2ba_pro, PBz5ba_e2ba, PBz5ba_e10ba, PBz8ba_e2ba, PBz8ba_e10ba, PBz10ba_ten, PBapc_e5ba,\
	PBapc_e8ba, PBe9starba, PBtenstar, PBz10ba_tenstar, PBPltsADP, PBPltse2, e2mtot, z2mtot, z5mtot, \
	e5mtot, z8mtot, e8mtot, z9mtot, e9mtot, e10mtot, z10mtot :"<<endl;
	for( int i = 0; i< solver->numV; i++)
	{
		cout<<domain.mesh[2][2].Cvals[i]<<", ";
	}

}
