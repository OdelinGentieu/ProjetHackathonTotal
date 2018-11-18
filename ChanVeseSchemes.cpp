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



double ChanVeseSchemes::Correction(const field& phi)
{
// Calcul de C1 et C2
double dom_plus=0., dom_moins=0, z_plus=0., z_moins=0., C1, C2;

	for (int i=0; i<phi.cols() ; i++)
	{
		for (int j=0; j<phi.rows() ; j++)
		{
			dom_plus += max (phi(i,j)/max(abs(phi(i,j)),1.E-16),0.);
			dom_moins -= max (-phi(i,j)/max(abs(phi(i,j)),1.E-16),0.);
			//le max au dénominateur sert à ne jamais diviser par 0
			//dom_plus est le nombre d'éléments positifs dans phi
			//dom_moins est le nombre d'éléments négatifs dans phi
			z_plus += max (_u0(i,j)*phi(i,j)/max(abs(phi(i,j)),1.E-16),0.);
			z_moins -= max (-_u0(i,j)*phi(i,j)/max(abs(phi(i,j)),1.E-16),0.);
			//z_plus est la valeur de l'intégrale de z sur tous les phi(i,j) positifs
			//z_moins est la valeur de l'intégrale de z sur tous les phi(i,j) négatifs
		}
	}
	C1 = z_plus/dom_plus;
	C2 = z_moins/dom_moins;
	if (C1 < C2) // On sait pas trop pourquoi mais c'était fait comme ça dans la version précédente, donc on a fait pareil
	{
		double temp = C2;
		C2 = C1; C1 = temp;
	}
	// Fin de Calcul de C1 et C2
	return (C1,C2);
}


field ChanVeseSchemes::ExplicitScheme(const field& Phi, const double dt,  const double mu, const double nu, \
	                                    const double l1, const double l2, const double C1, const double C2) const
{
	const double hx(1.), hy(1.);
	const double eta(1e-8);

	int nx_local(Phi.rows());
	int ny_local(Phi.cols());

	field phiinterior(nx_local,ny_local);

	double eps(3.);
	double diracij;

	// Application du Schema
	for (int i=1; i<nx_local; ++i)
	{
		for (int j=1; j<ny_local; ++j)
		{
			double firstterm  = (fdxplus(i,j,Phi,hx)*coeffA(i,j,Phi,hx,hy,eta) - fdxminus(i,j,Phi,hx)*coeffA(i-1,j,Phi,hx,hy,eta));
			double secondterm = (fdyplus(i,j,Phi,hy)*coeffB(i,j,Phi,hx,hy,eta) - fdyminus(i,j,Phi,hy)*coeffB(i,j-1,Phi,hx,hy,eta));
			double correc = -l1*(_u0(i-1,j-1)-C1)*(_u0(i-1,j-1)-C1) + l2*(_u0(i-1,j-1)-C2)*(_u0(i-1,j-1)-C2);
			diracij = eps/(Phi(i,j)*Phi(i,j)+eps*eps);
			phiinterior(i,j) = Phi(i,j) + dt*diracij*(mu*(firstterm+secondterm)- nu + correc);
		}
	}

	return phiinterior;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// // Valeur moyenne du domaine
// double ChanVeseSchemes::ComputeMeanValueOnDomain(const field& phi) const
// {
// 	field P = (phi>= 0).cast<double>();
// 	double size_domain = fmax(P.sum(),1.e-6);
// 	double zP = (_u0*P).sum();
// 	return 1./size_domain*zP;
// }
// // Valeur moyenne sur le domaine complémentaire
// double ChanVeseSchemes::ComputeMeanValueOnComplementaryDomain(const field& phi) const
// {
// 	field Q = (phi<= 0).cast<double>();
// 	double size_domain = fmax(Q.sum(),1.e-6);
// 	double zQ = (_u0*Q).sum();
// 	return 1./size_domain*zQ;
// }
// field ChanVeseSchemes::Correction(const field& phi, const double lambda1, const double lambda2) const
// {
// 	double Cmin(ComputeMeanValueOnDomain(phi)), Cmax(ComputeMeanValueOnComplementaryDomain(phi));
// 	if (Cmax < Cmin)
// 	{
// 		double temp = Cmin;
// 		Cmin = Cmax; Cmax = temp;
// 	}
// 	field correc_term_1 = (_u0-Cmin)*(_u0-Cmin);
// 	field correc_term_2 = (_u0-Cmax)*(_u0-Cmax);
// 	field correc_term = (-lambda1*correc_term_1+lambda2*correc_term_2);
// 	return correc_term;
// }


#endif
