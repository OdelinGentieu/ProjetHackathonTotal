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

	field dxplus  = ( CSXPshift(phi) - phi ) / (hx);
	field dxminus = ( phi - CSXMshift(phi) ) / (hx);
	field dyplus  = ( CSYPshift(phi) - phi ) / (hy);
	field dyminus = ( phi - CSYMshift(phi) ) / (hy);
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
	field correc_term = (-lambda1*correc_term_1+lambda2*correc_term_2);
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
	field transfophi(nx,ny);

	cout << "Taille nx et ny " << nx << " " << ny << endl;

	for (int i=2; i<nx-1; ++i)
	{
		for (int j=2; j<ny-1; ++j)
		{
			double firstAterm = fdxplus(i,j,phi,hx)/(hx * sqrt(pow(eta,2) + pow(fdxplus(i,j,phi, hx),2) + pow((fdycentral(i,j,phi, hy)/2.),2)));
			double secondAterm = fdxminus(i,j,phi,hx)/(hx * sqrt(pow(eta,2) + pow(fdxplus(i-1,j,phi, hx),2) + pow((fdycentral(i-1,j,phi, hy)/2.),2)));

			double firstBterm = fdyplus(i,j,phi, hy)/(hy * sqrt(pow(eta,2) + pow(fdyplus(i,j,phi, hy),2) + pow(fdxcentral(i,j,phi, hx),2)));
			double secondBterm = fdyminus(i,j,phi, hy)/(hx * sqrt(pow(eta,2) + pow(fdyplus(i,j-1,phi, hy),2) + pow(fdxcentral(i,j-1,phi, hx),2)));

			double eps(3.);
			double diracij;
			diracij=eps/(phi(i,j)*phi(i,j)+eps*eps);

			transfophi(i,j) = phi(i,j) + dt*diracij*mu*(firstAterm+secondAterm+firstBterm+secondBterm) + dt*diracij*(nu + correction(i,j));
		}
	}
  //
	// // Différences à gauche, à droite et centrées
	// field dxplus  = ( CSXPshift(phi) - phi ) / (hx);
	// field dxminus = ( phi - CSXMshift(phi) ) / (hx);
	// field dyplus  = ( CSYPshift(phi) - phi ) / (hy);
	// field dyminus = ( phi - CSYMshift(phi) ) / (hy);
	// field dxcentral = (dxplus+dxminus) / 2.;
	// field dycentral = (dyplus+dyminus) / 2.;
  //
	// field firstterm = dxplus/sqrt(eta*eta + dxplus*dxplus + dycentral*dycentral);
	// field secondterm = dyplus/sqrt(eta*eta + dyplus*dyplus + dxcentral*dxcentral);
  //
	// // Courbure
	// field curvature =  (firstterm-CSXMshift(firstterm))/hx + (secondterm-CSYMshift(secondterm))/hy ;
  //
	// // Dirac
	// field dirac = Dirac(phi);
  //
	// // Terme correctif
	// field correction = Correction(phi,l1,l2);
  //
	// return phi + dt*dirac*(mu*curvature-nu+correction);

	cout << "test 6" << endl;

	return transfophi;
}



#endif
