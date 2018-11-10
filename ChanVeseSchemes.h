// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _CHANVESESCHEMES_H

#include "Dense"
#include "Image.h"
#include "Util.h"
#include "LevelSet.h"

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class ChanVeseSchemes
{


private:
	// Image à segmenter
	field _u0;
	//double hx, hy;//, hy;



public:
	ChanVeseSchemes(Image* image);

	// Dérivées partielles
	field CSXPshift(const field& phi) const;
	field CSXMshift(const field& phi) const;
	field CSYPshift(const field& phi) const;
	field CSYMshift(const field& phi) const;
	field CSXMYPshift(const field& phi) const;
	field CSXPYMshift(const field& phi) const;
	field CSXMYMshift(const field& phi) const;

	// Discretisation du Dirac
	field Dirac(const field& phi) const;

	// Valeur moyenne du domaine
	double ComputeMeanValueOnDomain(const field& phi) const;

	// Valeur moyenne sur le domaine complémentaire
	double ComputeMeanValueOnComplementaryDomain(const field& phi) const;

	// Correction
	field Correction(const field& phi, const double lambda1, const double lambda2) const;

	// |V phi|
	field AbsGradPhi(const field& phi) const;

	// Schéma pour différences finis
	// Explicit Scheme
	field ExplicitScheme(const field& phi, const double dt,  const double mu, const double nu, const double l1, const double l2) const;

	inline double fdxplus(int i,int j,const field& phi, double hx) const
	{
		std::cout << "test 1 ter" <<std::endl;
		 return (phi(i+1,j)-phi(i,j))/hx;
	};

	inline double fdxminus(int i,int j, const field& phi, double hx) const
	{
		return (phi(i,j)-phi(i-1,j))/hx;
	};

	inline double fdyplus(int i,int j, const field& phi, double hy) const
	{
		(phi(i,j+1)-phi(i,j))/hy;
	};

	inline double fdyminus(int i,int j, const field& phi, double hy) const
	{
		return (phi(i,j)-phi(i,j-1))/hy;
	};

	inline double fdxcentral(int i,int j, const field& phi, double hx) const
	{
		return (fdxplus(i,j,phi, hx)+fdxminus(i,j, phi,hx)) / 2.;
	};

	inline double fdycentral(int i,int j,const field& phi,double hy) const
	{
		return (fdyplus(i,j,phi, hy)+fdyminus(i,j,phi, hy)) / 2.;
	};
};

#define _CHANVESESCHEMES_H
#endif
