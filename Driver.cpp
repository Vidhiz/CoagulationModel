#include <iostream>
#include <fstream>
#include <cmath>
#include "ODEsolver.h"
#include "Domain.h"

using namespace std;
using namespace Reactions;

// signature of function to display output concentrations:
void displayC(ODESolver *, Domain&, double );

// function that myrounds to desired precision
double myround(double f,int prec);

/////////////////////////////
int main(int argc, char *argv[])
{
	ODESolver *solver = new ODESolver(1.0000000000000002e-002,1.0000000000000000E-003);

	Domain domain(80,20,128,32);
	Mesh m = domain.get_mesh();

	// open files and load m with data 
	ifstream myfilepba ("/home/vidhizala/Documents/karincode/Data/befor.Pba.0032.0001.dat", ios::in);
	ifstream myfilepsea ("/home/vidhizala/Documents/karincode/Data/befor.Psea.0032.0001.dat",ios::in);
	ifstream myfilez2 ("/home/vidhizala/Documents/karincode/Data/befor.z2ba.0032.0001.dat ",ios::in);
	ifstream myfilee2 ("/home/vidhizala/Documents/karincode/Data/befor.e2ba.0032.0001.dat",ios::in);
	ifstream myfilez5 ("/home/vidhizala/Documents/karincode/Data/befor.z5ba.0032.0001.dat",ios::in);
	ifstream myfilee5 ("/home/vidhizala/Documents/karincode/Data/befor.e5ba.0032.0001.dat",ios::in);
	ifstream myfilez8 ("/home/vidhizala/Documents/karincode/Data/befor.z8ba.0032.0001.dat",ios::in);
	ifstream myfilee8 ("/home/vidhizala/Documents/karincode/Data/befor.e8ba.0032.0001.dat",ios::in);
	ifstream myfilez9 ("/home/vidhizala/Documents/karincode/Data/befor.z9ba.0032.0001.dat",ios::in);
	ifstream myfilee9 ("/home/vidhizala/Documents/karincode/Data/befor.e9ba.0032.0001.dat",ios::in);
	ifstream myfilez10 ("/home/vidhizala/Documents/karincode/Data/befor.z10ba.0032.0001.dat",ios::in);
	ifstream myfilee10 ("/home/vidhizala/Documents/karincode/Data/befor.e10ba.0032.0001.dat",ios::in);
	ifstream myfileten ("/home/vidhizala/Documents/karincode/Data/befor.ten.0032.0001.dat",ios::in);
	ifstream myfilepro ("/home/vidhizala/Documents/karincode/Data/befor.pro.0032.0001.dat",ios::in);
	ifstream myfilez2bapro ("/home/vidhizala/Documents/karincode/Data/befor.z2ba_pro.0032.0001.dat ",ios::in);
	ifstream myfilez5bae2ba ("/home/vidhizala/Documents/karincode/Data/befor.z5ba_e2ba.0032.0001.dat",ios::in);
	ifstream myfilez5bae10ba ("/home/vidhizala/Documents/karincode/Data/befor.z5ba_e10ba.0032.0001.dat",ios::in);
	ifstream myfilez8bae2ba ("/home/vidhizala/Documents/karincode/Data/befor.z8ba_e2ba.0032.0001.dat",ios::in);
	ifstream myfilez8bae10ba ("/home/vidhizala/Documents/karincode/Data/befor.z8ba_e10ba.0032.0001.dat",ios::in);
	ifstream myfilez10baten ("/home/vidhizala/Documents/karincode/Data/befor.z10ba_ten.0032.0001.dat",ios::in);
	ifstream myfileapce5ba ("/home/vidhizala/Documents/karincode/Data/befor.apc_e5ba.0032.0001.dat",ios::in);
	ifstream myfileapce8ba ("/home/vidhizala/Documents/karincode/Data/befor.apc_e8ba.0032.0001.dat",ios::in);
	ifstream myfilee9starba ("/home/vidhizala/Documents/karincode/Data/befor.e9starba.0032.0001.dat",ios::in);
	ifstream myfilepbtenstar ("/home/vidhizala/Documents/karincode/Data/befor.tenstar.0032.0001.dat",ios::in);
	ifstream myfilez10batenstar ("/home/vidhizala/Documents/karincode/Data/befor.z10ba_tenstar.0032.0001.dat",ios::in);

	int count;
	double tmp;
	myfilepba >> count;
	// for each point in the domain :
	for(int j = 0; j<domain.ny; j++)
	{
		for(int k = 0; k<domain.nx; k++)
		{
			
	        myfilepba >> tmp;
	        m[j][k].Pba = tmp;
	        //cout << m[j][k].Pba<<"\n";

	        myfilepsea >> tmp;
	        m[j][k].Psea = tmp;

	        myfilez2 >> tmp;
	        m[j][k].Cvals[PBz2ba] = tmp;

	        myfilee2 >> tmp;
	        m[j][k].Cvals[PBe2ba] = tmp;

	        myfilez5 >> tmp;
	        m[j][k].Cvals[PBz5ba] = tmp;

	        myfilee2 >> tmp;
	        m[j][k].Cvals[PBe5ba] = tmp;

	        myfilez8 >> tmp;
	        m[j][k].Cvals[PBz8ba] = tmp;

	        myfilee8 >> tmp;
	        m[j][k].Cvals[PBe8ba] = tmp;

	        myfilez9 >> tmp;
	        m[j][k].Cvals[PBz9ba] = tmp;

	        myfilee9 >> tmp;
	        m[j][k].Cvals[PBe9ba] = tmp;

	        myfilez10 >> tmp;
	        m[j][k].Cvals[PBz10ba] = tmp;

	        myfilee10 >> tmp;
	        m[j][k].Cvals[PBe10ba] = tmp;

	        myfileten >> tmp;
	        m[j][k].Cvals[PBten] = tmp;

	        myfilepro >> tmp;
	        m[j][k].Cvals[PBpro] = tmp;

	        myfilez2bapro >> tmp;
	        m[j][k].Cvals[PBz2ba_pro] = tmp;

	        myfilez5bae2ba >> tmp;
	        m[j][k].Cvals[PBz5ba_e2ba] = tmp;

	        myfilez5bae10ba >> tmp;
	        m[j][k].Cvals[PBz5ba_e10ba] = tmp;

	        myfilez8bae2ba >> tmp;
	        m[j][k].Cvals[PBz8ba_e2ba] = tmp;
	        
	        myfilez8bae10ba >> tmp;
	        m[j][k].Cvals[PBz8ba_e10ba] = tmp;
	        
	        myfilez10baten >> tmp;
	        m[j][k].Cvals[PBz10ba_ten] = tmp;
	        
	        myfileapce5ba >> tmp;
	        m[j][k].Cvals[PBapc_e5ba] = tmp;
	        
	        myfileapce8ba >> tmp;
	        m[j][k].Cvals[PBapc_e8ba] = tmp;

	        myfilee9starba>> tmp;
	        m[j][k].Cvals[PBe9starba] = tmp;

	        myfilepbtenstar >> tmp;
	        m[j][k].Cvals[PBtenstar] = tmp;

	        myfilez10batenstar >> tmp;
	        m[j][k].Cvals[PBz10ba_tenstar] = tmp;
		   
		}
	}

	double current_time = 0; 

	// begin time tracking:
	time_t start = time(NULL);

	current_time = solver->RK2solve(domain);

	cout<<"time elapsed="<<time(NULL) - start<<endl;;

	displayC(solver, domain, current_time);


	return 0;
}

void displayC(ODESolver  *solver, Domain &domain, double current_time)
{
	/*cout<<"\nComputation complete. Here are the concentrations for chemicals at a point (x,y) in time "
	<<current_time<<" s"<<endl;
	cout<<" PBz2ba, PBe2ba, PBz5ba, PBe5ba, PBz8ba, PBe8ba, PBz9ba, PBe9ba,PBz10ba, PBe10ba, PBten,\
	PBpro, PBz2ba_pro, PBz5ba_e2ba, PBz5ba_e10ba, PBz8ba_e2ba, PBz8ba_e10ba, PBz10ba_ten, PBapc_e5ba,\
	PBapc_e8ba, PBe9starba, PBtenstar, PBz10ba_tenstar, PBPltsADP, PBPltse2, e2mtot, z2mtot, z5mtot, \
	e5mtot, z8mtot, e8mtot, z9mtot, e9mtot, e10mtot, z10mtot :"<<endl;
	for( int i = 0; i< solver->numV; i++)
	{
		cout<<domain.mesh[2][2].Cvals[i]<<", ";
	}*/

	ifstream myfilez2 ("/home/vidhizala/Documents/karincode/Data/after.z2ba.0032.0002.dat ",ios::in);
	ifstream myfilee2 ("/home/vidhizala/Documents/karincode/Data/after.e2ba.0032.0002.dat",ios::in);
	ifstream myfilez5 ("/home/vidhizala/Documents/karincode/Data/after.z5ba.0032.0002.dat",ios::in);
	ifstream myfilee5 ("/home/vidhizala/Documents/karincode/Data/after.e5ba.0032.0002.dat",ios::in);
	ifstream myfilez8 ("/home/vidhizala/Documents/karincode/Data/after.z8ba.0032.0002.dat",ios::in);
	ifstream myfilee8 ("/home/vidhizala/Documents/karincode/Data/after.e8ba.0032.0002.dat",ios::in);
	ifstream myfilez9 ("/home/vidhizala/Documents/karincode/Data/after.z9ba.0032.0002.dat",ios::in);
	ifstream myfilee9 ("/home/vidhizala/Documents/karincode/Data/after.e9ba.0032.0002.dat",ios::in);
	ifstream myfilez10 ("/home/vidhizala/Documents/karincode/Data/after.z10ba.0032.0002.dat",ios::in);
	ifstream myfilee10 ("/home/vidhizala/Documents/karincode/Data/after.e10ba.0032.0002.dat",ios::in);
	ifstream myfileten ("/home/vidhizala/Documents/karincode/Data/after.ten.0032.0002.dat",ios::in);
	ifstream myfilepro ("/home/vidhizala/Documents/karincode/Data/after.pro.0032.0002.dat",ios::in);
	ifstream myfilez2bapro ("/home/vidhizala/Documents/karincode/Data/after.z2ba_pro.0032.0002.dat ",ios::in);
	ifstream myfilez5bae2ba ("/home/vidhizala/Documents/karincode/Data/after.z5ba_e2ba.0032.0002.dat",ios::in);
	ifstream myfilez5bae10ba ("/home/vidhizala/Documents/karincode/Data/after.z5ba_e10ba.0032.0002.dat",ios::in);
	ifstream myfilez8bae2ba ("/home/vidhizala/Documents/karincode/Data/after.z8ba_e2ba.0032.0002.dat",ios::in);
	ifstream myfilez8bae10ba ("/home/vidhizala/Documents/karincode/Data/after.z8ba_e10ba.0032.0002.dat",ios::in);
	ifstream myfilez10baten ("/home/vidhizala/Documents/karincode/Data/after.z10ba_ten.0032.0002.dat",ios::in);
	ifstream myfileapce5ba ("/home/vidhizala/Documents/karincode/Data/after.apc_e5ba.0032.0002.dat",ios::in);
	ifstream myfileapce8ba ("/home/vidhizala/Documents/karincode/Data/after.apc_e8ba.0032.0002.dat",ios::in);
	ifstream myfilee9starba ("/home/vidhizala/Documents/karincode/Data/after.e9starba.0032.0002.dat",ios::in);
	ifstream myfilepbtenstar ("/home/vidhizala/Documents/karincode/Data/after.tenstar.0032.0002.dat",ios::in);
	ifstream myfilez10batenstar ("/home/vidhizala/Documents/karincode/Data/after.z10ba_tenstar.0032.0002.dat",ios::in);

	int count, flag = 0;
	bool loop = true;
	double tmp;
	Mesh m = domain.get_mesh();
	myfilez2 >> count;
	// for each point in the domain :
	for(int j = 0; j<domain.ny && loop; j++)
	{
		for(int k = 0; k<domain.nx && loop; k++)
		{
			
	        myfilez2 >> tmp;
	        if(myround(m[j][k].Cvals[PBz2ba],4) != myround(tmp,4))
 			{
	        	cout<<"\n3\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }

	        myfilee2 >> tmp;
	        if(myround(m[j][k].Cvals[PBe2ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n4\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }

	        myfilez5 >> tmp;
	        if(myround(m[j][k].Cvals[PBz5ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n5\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }

	        myfilee2 >> tmp;
	        if(myround(m[j][k].Cvals[PBe5ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n4\nbreaking! ";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez8 >> tmp;
	        if(myround(m[j][k].Cvals[PBz8ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n7\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilee8 >> tmp;
	        if(myround(m[j][k].Cvals[PBe8ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n8\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez9 >> tmp;
	        if(myround(m[j][k].Cvals[PBz9ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n9\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilee9 >> tmp;
	        if(myround(m[j][k].Cvals[PBe9ba] ,4)!= myround(tmp,4))
	         {
	        	cout<<"\n10\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez10 >> tmp;
	        if(myround(m[j][k].Cvals[PBz10ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n11\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilee10 >> tmp;
	        if(myround(m[j][k].Cvals[PBe10ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n12\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	       
	        myfileten >> tmp;
	        if(myround(m[j][k].Cvals[PBten],4) != myround(tmp,4))
	         {
	        	cout<<"\n13\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilepro >> tmp;
	        if(myround(m[j][k].Cvals[PBpro],4) != myround(tmp,4))
	         {
	        	cout<<"\n14\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez2bapro >> tmp;
	        if(myround(m[j][k].Cvals[PBz2ba_pro] ,4)!= myround(tmp,4))
	         {
	        	cout<<"\n15\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez5bae2ba >> tmp;
	        if(myround(m[j][k].Cvals[PBz5ba_e2ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n16\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez5bae10ba >> tmp;
	        if(myround(m[j][k].Cvals[PBz5ba_e10ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n17\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez8bae2ba >> tmp;
	        if(myround(m[j][k].Cvals[PBz8ba_e2ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n18\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez8bae10ba >> tmp;
	        if(myround(m[j][k].Cvals[PBz8ba_e10ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n19\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez10baten >> tmp;
	        if(myround(m[j][k].Cvals[PBz10ba_ten],4) != myround(tmp,4))
	         {
	        	cout<<"\n20\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfileapce5ba >> tmp;
	        if(myround(m[j][k].Cvals[PBapc_e5ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n21\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfileapce8ba >> tmp;
	        if(myround(m[j][k].Cvals[PBapc_e8ba],4) != myround(tmp,4))
	         {
	        	cout<<"\n22\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilee9starba>> tmp;
	        if(myround(m[j][k].Cvals[PBe9starba],4) != myround(tmp,4))
	         {
	        	cout<<"\n23\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilepbtenstar >> tmp;
	        if(myround(m[j][k].Cvals[PBtenstar],4) != myround(tmp,4))
	         {
	        	cout<<"\n24\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
	        myfilez10batenstar >> tmp;
	        if(myround(m[j][k].Cvals[PBz10ba_tenstar],4) != myround(tmp,4))
		   	 {
	        	cout<<"\n25\nbreaking!";
	        	flag = 1; loop = false;
	        	break;
	        }
		}
	}

	if(flag == 0)
	{
		cout<<"\n Verified.";
	}
	cout<<"\n Test complete!";

}

double myround(double f,int prec)
{
    return (double) (floor(f*(1.0f/prec) + 0.5)/(1.0f/prec));
}
