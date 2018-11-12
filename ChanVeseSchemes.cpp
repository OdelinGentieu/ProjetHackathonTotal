// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18



#ifndef _CHANVESESCHEMES_CPP
#define _CHANVESESCHEMES_CPP


#include "ChanVeseSchemes.h"
#include "float.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


ChanVeseSchemes::ChanVeseSchemes (Image* image) : _u0(image->GetImage()) {}

// Dérivées partielles
field ChanVeseSchemes::CSXPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXP(nx,ny); CSXP.leftCols(ny-1) = phi.rightCols(ny-1); CSXP.col(ny-1) = phi.col(0);
	return CSXP;
}
field ChanVeseSchemes::CSXMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXM(nx,ny); CSXM.rightCols(ny-1) = phi.leftCols(ny-1); CSXM.col(0) = phi.col(ny-1);
	return CSXM;
}
field ChanVeseSchemes::CSYPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSYP(nx,ny); CSYP.topRows(nx-1) = phi.bottomRows(nx-1); CSYP.row(nx-1) = phi.row(0);
	return CSYP;
}
field ChanVeseSchemes::CSYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSYM(nx,ny); CSYM.bottomRows(nx-1) = phi.topRows(nx-1); CSYM.row(0) = phi.row(nx-1);
	return CSYM;
}

/////////////////////////////////////////////////////////////////////////////////////////:///

field ChanVeseSchemes::CSXMYPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXMYP(nx,ny); CSXMYP = CSYPshift(phi); CSXMYP = CSXMshift(CSXMYP);
	return CSXMYP;
}
field ChanVeseSchemes::CSXPYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXPYM(nx,ny); CSXPYM = CSYMshift(phi); CSXPYM = CSXPshift(CSXPYM);
	return CSXPYM;
}
field ChanVeseSchemes::CSXMYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXMYM(nx,ny); CSXMYM = CSYMshift(phi); CSXMYM = CSXMshift(CSXMYM);
	return CSXMYM;
}

// Discretisation du dirac
field ChanVeseSchemes::Dirac(const field& phi) const
{
	double eps(3.);
	return eps/(phi*phi+eps*eps);
}

// Discretisation de la valeur absolue du gradient
field ChanVeseSchemes::AbsGradPhi(const field& phi) const
{
	const double hx(1.), hy(1.0);

	field dxplus    = ( CSXPshift(phi) - phi ) / (hx);
	field dxminus   = ( phi - CSXMshift(phi) ) / (hx);
	field dyplus    = ( CSYPshift(phi) - phi ) / (hy);
	field dyminus   = ( phi - CSYMshift(phi) ) / (hy);
	field dxcentral = (dxplus+dxminus) / 2.;
	field dycentral = (dyplus+dyminus) / 2.;

	return sqrt(dxcentral*dxcentral + dycentral*dycentral);
}


// Valeur moyenne du domaine
double ChanVeseSchemes::ComputeMeanValueOnDomain(const field& phi) const
{
	field P = (phi>= 0).cast<double>();
	double size_domain = fmax(P.sum(),1.e-6);
	double zP = (_u0*P).sum();
	return zP/size_domain;
}
// Valeur moyenne sur le domaine complémentaire
double ChanVeseSchemes::ComputeMeanValueOnComplementaryDomain(const field& phi) const
{
	field Q = (phi<= 0).cast<double>();
	double size_domain = fmax(Q.sum(),1.e-6);
	double zQ = (_u0*Q).sum();
	return zQ/size_domain;
}

field ChanVeseSchemes::Correction(const field& phi, const double lambda1, const double lambda2) const
{
	double Cmin(ComputeMeanValueOnDomain(phi)), Cmax(ComputeMeanValueOnComplementaryDomain(phi));
	if (Cmax < Cmin)
	{
		double temp = Cmin;
		Cmin = Cmax; Cmax = temp;
	}
	field correc_term_1 = (_u0-Cmin)*(_u0-Cmin);
	field correc_term_2 = (_u0-Cmax)*(_u0-Cmax);
	field correc_term   = (-lambda1*correc_term_1+lambda2*correc_term_2);
	return correc_term;
}

field ChanVeseSchemes::ExplicitScheme(const field& phi, const double dt,  const double mu, const double nu, const double l1, const double l2) const
{

	double hx(1.0);
	double hy(1.0);
	const double eta(1e-8);

	int nx(phi.rows());
	int ny(phi.cols());

	field correction = Correction(phi,l1,l2);

	cout << "Taille nx et ny " << nx << " " << ny << endl;


	// DEBUT GROS TABLEAU

	field GrosPhi(nx+2,ny+2);
	field newphi(nx,ny);

	cout << "Taille new " << newphi.rows() << " " << newphi.cols() << " et l'enrobé " << GrosPhi.rows() << " " << GrosPhi.cols() << endl;
	cout << GrosPhi(0,0) << endl;

	for (int j=1; j<ny+1; ++j)
	{
		GrosPhi(0,j)    = phi(nx-1,j-1);
		GrosPhi(nx+1,j) = phi(0,j-1);
	}

	for (int i=1; i<nx+1; ++i)
	{
		GrosPhi(i,0)    = phi(i-1,ny-1);
		GrosPhi(i,ny+1) = phi(i-1,0);
		for (int j=1; j<ny+1; ++j)
		{
			GrosPhi(i,j)  = phi(i-1,j-1);
		}
	}

	cout << "test1" << endl;
	//FIN GrosPhi
	for (int i=1; i<nx+1; ++i)
	{
		for (int j=1; j<ny+1; ++j)
		{
			double firstterm  = fdxplus(i,j,GrosPhi,hx)*coeffA(i,j,GrosPhi,hx,hy,eta) + fdxminus(i,j,GrosPhi,hx)*coeffA(i-1,j,GrosPhi,hx,hy,eta);
			double secondterm = fdyplus(i,j,GrosPhi,hy)*coeffB(i,j,GrosPhi,hx,hy,eta) + fdyminus(i,j,GrosPhi,hy)*coeffB(i,j-1,GrosPhi,hx,hy,eta);
			cout << "test2" << endl;
			double eps(3.);
			double diracij;
			diracij = eps/(GrosPhi(i,j)*GrosPhi(i,j)+eps*eps);
			cout << "test3" << endl;
			newphi(i-1,j-1) = GrosPhi(i,j) + dt*diracij*mu*(firstterm+secondterm) + dt*diracij*(nu + correction(i-1,j-1));

			cout << "test4" << endl;
		}
	}



	return newphi;
}



#endif
