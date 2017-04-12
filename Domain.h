#ifndef _REACTIONS_DOMAIN_H_
#define _REACTIONS_DOMAIN_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "ODEsolver.h"
#include "Point.h"

namespace Reactions{

typedef std::vector<std::vector<Point>> Mesh;

/// Domain
class Domain{
	public:
		/// Default Constructor
		Domain()
		{
			std::cout<<"Created Domain object.\n";
			xlen = 2.0;
			ylen = 1.0;
			nx = 40;
			ny = 20;

			// calculate dx and dy based on user inputs
			dx = xlen/nx;
			dy = ylen/ny;

			// Set up mesh size. (HEIGHT x WIDTH)
			mesh.resize(ny);
			for (int i = 0; i < ny; ++i)
				mesh[i].resize(nx);
	
		}

		/// Parametrized constructor
		Domain(double xlen, double ylen, double nx, double ny)
		{
			Domain::xlen = xlen;
			Domain::ylen = ylen;
			Domain::nx = nx;
			Domain::ny = ny;

			// Set up mesh size. (HEIGHT x WIDTH)
			mesh.resize(ny);
			for (int i = 0; i < ny; ++i)
				mesh[i].resize(nx);
			

		}

		/// Destructor
		~Domain()
		{
			std::cout<<"Destroyed Domain object.\n";
		}	

		// mesh
		Mesh mesh;

		/// x - length
		double xlen;

		/// y - length
		double ylen;

		/// points in x
		double nx;

		/// points in y
		double ny;

		// grid spacing in x
		double dx;

		// grid spacing in y
		double dy;

		/// get variables:
		double& get_nx()
		{
			return this->nx;
		}
		double& get_ny()
		{
			return this->ny;
		}
		double& get_xlen()
		{
			return this->xlen;
		}
		double& get_ylen()
		{
			return this->ylen;
		}
		double& get_dx()
		{
			return this->dx;
		}
		double& get_dy()
		{
			return this->dy;
		}
		Mesh& get_mesh()
		{
			return this->mesh;
		}

		/// set variables:
		void set_nx(double n)
		{
			this->nx = n;
		}
		void set_ny(double n)
		{
			this->ny = n;
		}
		void set_xlen(double len)
		{
		    this->xlen = len;
		}
		void set_ylen(double len)
		{
			this->ylen = len;
		}
		void set_dx(double d)
		{
			this->dx = d;
		}
		void set_dy(double d)
		{
			this->dy = d;
		}
		void set_mesh(Mesh &m)
		{
			this->mesh = m;
		}

	   void updateTotal()
	   {
	   		for(int j = 0; j<ny; j++)
				for(int k = 0; k<nx;k++)
	   			{
	   				mesh[j][k].Cvals[e2mtot] = mesh[j][k].Cvals[PBe2ba] + mesh[j][k].Cvals[PBz5ba_e2ba] + mesh[j][k].Cvals[PBz8ba_e2ba]; 
	   				mesh[j][k].Cvals[z2mtot] = mesh[j][k].Cvals[PBz2ba] + mesh[j][k].Cvals[PBz2ba_pro];
	   				mesh[j][k].Cvals[ze5mtot] = mesh[j][k].Cvals[PBz5ba] + mesh[j][k].Cvals[PBz5ba_e10ba] +  mesh[j][k].Cvals[PBz5ba_e2ba] + mesh[j][k].Cvals[PBe5ba] + mesh[j][k].Cvals[PBapc_e5ba] + mesh[j][k].Cvals[PBpro] + mesh[j][k].Cvals[PBz2ba_pro];
	   				mesh[j][k].Cvals[ze8mtot] = mesh[j][k].Cvals[PBz8ba] + mesh[j][k].Cvals[PBz8ba_e2ba] + mesh[j][k].Cvals[PBz8ba_e10ba] + mesh[j][k].Cvals[PBe8ba] + mesh[j][k].Cvals[PBapc_e8ba] + mesh[j][k].Cvals[PBten] + mesh[j][k].Cvals[PBz10ba_ten] + mesh[j][k].Cvals[PBtenstar] + mesh[j][k].Cvals[PBz10ba_tenstar];
	   				mesh[j][k].Cvals[ze9mtot] = mesh[j][k].Cvals[PBz9ba] + mesh[j][k].Cvals[PBe9ba] + mesh[j][k].Cvals[PBz10ba_ten] + mesh[j][k].Cvals[PBten];
	   				mesh[j][k].Cvals[ze10mtot] = mesh[j][k].Cvals[PBz10ba] + mesh[j][k].Cvals[PBz10ba_tenstar] + mesh[j][k].Cvals[PBe10ba] + mesh[j][k].Cvals[PBpro] + mesh[j][k].Cvals[PBz8ba_e10ba] + mesh[j][k].Cvals[PBz5ba_e10ba] + mesh[j][k].Cvals[PBz10ba_ten] + mesh[j][k].Cvals[PBz2ba_pro];

	   			}
	   }
  
		
};

} // namespace Reactions
#endif //_REACTIONS_DOMAIN_H_