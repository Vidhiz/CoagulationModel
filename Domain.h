
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Point.h"


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
			
			// fill mesh with coordinates


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

		/// get and set variables:
		double get_nx()
		{
			return nx;
		}
		double get_ny()
		{
			return ny;
		}
		double get_xlen()
		{
			return xlen;
		}
		double get_ylen()
		{
			return ylen;
		}
		double get_dx()
		{
			return dx;
		}
		double get_dy()
		{
			return dy;
		}
		Mesh get_mesh()
		{
			return mesh;
		}
		
};