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

ChanVeseSchemes::ChanVeseSchemes (Image* image) : _u0(image->GetImage())
{
	//création de _u0_v
	_u0_v.resize(_u0.rows());

	for (int i=0 ; i< _u0.rows(); i++)
	{
		_u0_v[i].resize(_u0.cols());
		for (int j=0 ; j< _u0.cols(); j++)
		{
			_u0_v[i][j]=_u0(i,j);
		}
	}
}


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
	// diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()/(phi.rows()*phi.cols());
	double d=0.;

	for (int i=0; i<newphi_v.size(); i++)
	{
		for (int j=0; j<newphi_v[0].size(); j++)
		{
			// if (phi_v[i][j] >= 0. && newphi_v[i][j] <0.)
			// {
			// 	d+= 4;
			// };
			// if (phi_v[i][j]< 0. && newphi_v[i][j] >=0.)
			// {
			// 	d+=4;
			// };

			d+= 4.*max( -phi_v[i][j]*newphi_v[i][j]/max(-phi_v[i][j]*newphi_v[i][j], 1.E-16) , 0. );
			//Petite modif par rapport au code d'origine : lorsque phi ou newphi était nul et l'autre étais négatif, d était incrémenté de 4.
		}
	}
	double diff = sqrt(d)/(newphi_v.size()*newphi_v[0].size());



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
	const double hx(1.), hy(1.0);
	const double eta(1e-8);

	int nx(phi_v.size());
	int ny(phi_v[0].size());

	std::vector< std::vector<double>> GrosPhi_v;
	GrosPhi_v.resize(nx+2);
	for (int i=0;i< nx+2;i++) { GrosPhi_v[i].resize(ny+2); }

	std::vector< std::vector<double>> newphi_v;
	newphi_v.resize(nx);
	for (int i=0;i<nx ;i++) { newphi_v[i].resize(ny); }


	// Calcul de C1 et C2
	double dom_plus=0., dom_moins=0, z_plus=0., z_moins=0., C1, C2;

	for (int i=0; i<nx ; i++)
	{
		for (int j=0; j<ny ; j++)
		{
			dom_plus += max(phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			dom_moins-= max(-phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			//le max au dénominateur sert à ne jamais diviser par 0
			//dom_plus est le nombre d'éléments positifs dans phi
			//dom_moins est le nombre d'éléments négatifs dans phi
			z_plus   += max(_u0_v[i][j]*phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			z_moins  -= max(-_u0_v[i][j]*phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			//z_plus est la valeur de l'intégrale de z sur tous les phi(i,j) positifs
			//z_moins est la valeur de l'intégrale de z sur tous les phi(i,j) négatifs
		}
	}
	C1 = z_plus/dom_plus;
	C2 = z_moins/dom_moins;

	// On sait pas trop pourquoi mais c'était fait comme ça dans la version précédente, donc on a fait pareil
	double temp = C2;
	C2 = max(C1,C2);
	C1 = min(C1,temp);

	//DEBUT Grosphi_v

	for (int j=1; j<ny+1; ++j)
	{
		GrosPhi_v[0][j]    = phi_v[0][j-1];
		GrosPhi_v[nx+1][j] = phi_v[nx-1][j-1];
	}

	for (int i=1; i<nx+1; ++i)
	{
		GrosPhi_v[i][0]    = phi_v[i-1][0];
		GrosPhi_v[i][ny+1] = phi_v[i-1][ny-1];
		for (int j=1; j<ny+1; ++j)
		{
			GrosPhi_v[i][j]  = phi_v[i-1][j-1];
		}
	}

	GrosPhi_v[0][0]        = phi_v[0][0];
	GrosPhi_v[0][ny+1]     = phi_v[0][ny-1];
	GrosPhi_v[nx+1][0]     = phi_v[nx-1][0];
	GrosPhi_v[nx+1][ny+1]  = phi_v[nx-1][ny-1];

	//FIN GrosPhi_v

	for (int i=1; i<nx+1; ++i)
	{
		for (int j=1; j<ny+1; ++j)
		{
			double firstterm   = (fdxplus(i,j,GrosPhi_v,hx)*coeffA(i,j,GrosPhi_v,hx,hy,eta) - fdxminus(i,j,GrosPhi_v,hx)*coeffA(i-1,j,GrosPhi_v,hx,hy,eta))/hx;
			double secondterm  = (fdyplus(i,j,GrosPhi_v,hy)*coeffB(i,j,GrosPhi_v,hx,hy,eta) - fdyminus(i,j,GrosPhi_v,hy)*coeffB(i,j-1,GrosPhi_v,hx,hy,eta))/hy;
			double eps(3.);
			double diracij;
			double correc      = -l1*(_u0(i-1,j-1)-C1)*(_u0(i-1,j-1)-C1) + l2*(_u0(i-1,j-1)-C2)*(_u0(i-1,j-1)-C2);
			diracij            = eps/(GrosPhi_v[i][j]*GrosPhi_v[i][j]+eps*eps);
			newphi_v[i-1][j-1] = GrosPhi_v[i][j] + dt*diracij*(mu*(firstterm+secondterm)- nu + correc);
		}
	}

	return newphi_v;

	////////////////// OLD STUFF ///////////////////

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
	//
	// newphi = phi + dt*dirac*(mu*curvature-nu+correction);
	//
	// cout << newphi(0,0) << " " << newphi(0,1) << " " << newphi(1,0) << " " << newphi(1,1) << endl;
	// cout << "Dirac " << dirac(0,0) << endl;
	//
	// return newphi;


}


#endif
