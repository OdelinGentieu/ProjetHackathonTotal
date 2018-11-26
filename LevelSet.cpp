#ifndef LEVELSET_CPP
#define LEVELSET_CPP

#include "LevelSet.h"

/*
 * redistancing max_iter times
 *
 */

void LevelSet::redistancing(const int max_iter)
{
  std::vector<std::vector<double>> phi0_v = _phi_v;
  std::vector<std::vector<double>> phin_v = _phi_v;

    int iter = 0, pourcentage =0;

    while (iter<max_iter)
      {
        firstOrderUpwind_v(_phi_v, phi0_v, .5*1);
        phin_v = _phi_v;

      	if ((floor(iter*100/max_iter)>=pourcentage)&&(max_iter>20))
      	  {
      	    pourcentage=pourcentage+1;

      	    int i_barre;
      	    printf( "[" );
      	    for (i_barre = 0 ; i_barre <= pourcentage ; i_barre += 2) printf( "*" );
      	    for (; i_barre <= 100 ; i_barre += 2 ) printf( "-" );
      	    printf( "] %3d %%", pourcentage );

      	    for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );

      	    fflush(stdout );

      	  }

        iter++;
      }

    if (max_iter > 20) {std::cout << std::endl;}
}

/*
 * redistancing until equilibrium
 *
 */
 void LevelSet::redistancing(const double eps)
 {
     field phi0 = _phi;
     field phin = _phi;

     int iter = 0; double err = 2*eps;
     while (err>eps)
     {
         firstOrderUpwind(_phi, phi0, .5*1);
         err = 2*(phin-_phi).abs().maxCoeff()/1;
         phin = _phi;
         iter++;
     }
 }

/*
 * first order upwind_v
 *
 */
void LevelSet::firstOrderUpwind_v(std::vector<std::vector<double>>& phi_v, const std::vector<std::vector<double>>& phi0_v, const double dt)
{
    std::vector<std::vector<double>> dxm = diff_m_v(phi_v, 0),
                                     dxp = diff_p_v(phi_v, 0),
                                     dym = diff_m_v(phi_v, 1),
                                     dyp = diff_p_v(phi_v, 1);

    std::vector<std::vector<double>> G = godunov_v(dxm, dxp, dym, dyp);
    std::vector<std::vector<double>> phi1_v = phi_v;

    for (int  i=0; i<phi_v.size(); ++i)
        for (int  j=0; j<phi_v[0].size(); ++j)
        {
            int mini(i+1); if (phi_v.size()-1 < mini) {mini = phi_v.size()-1;}
            int minj(j+1); if (phi_v[0].size()-1 < minj) {minj = phi_v[0].size()-1;}

            double  pij = phi0_v[i][j],
                    pimj = phi0_v[std::max(i-1,0)][j], pipj = phi0_v[mini][j],
                    pijm = phi0_v[i][std::max(j-1,0)], pijp = phi0_v[i][minj];
            double delta = boundary(pij, pimj, pipj, pijm, pijp);
            double Dij = .5*sgn(pij)*1; // Only for masks !

            phi1_v[i][j] = phi_v[i][j]-dt*delta*sgn(phi0_v[i][j])*G[i][j];
            phi1_v[i][j]-= dt*(1-delta)*( sgn(phi0_v[i][j])*fabs(phi_v[i][j])-Dij )/1;
        }
    phi_v = phi1_v;
}

/*
 * first order upwind
 *
 */
void LevelSet::firstOrderUpwind(field& phi, const field& phi0, const double dt)
{
    field dxm = diff_m(phi, 0),
          dxp = diff_p(phi, 0),
          dym = diff_m(phi, 1),
          dyp = diff_p(phi, 1);

    field G = godunov(dxm, dxp, dym, dyp);

    field phi1 = phi;
    for (int  i=0; i<phi.rows(); ++i)
        for (int  j=0; j<phi.cols(); ++j)
        {
            int mini(i+1); if (phi.rows()-1 < mini) {mini = phi.rows()-1;}
            int minj(j+1); if (phi.cols()-1 < minj) {minj = phi.cols()-1;}

            double  pij = phi0(i,j),
                    pimj = phi0(std::max(i-1,0),j), pipj = phi0(mini,j),
                    pijm = phi0(i,std::max(j-1,0)), pijp = phi0(i,minj);
            double delta = boundary(pij, pimj, pipj, pijm, pijp);
            double Dij = .5*sgn(pij)*1; // Only for masks !

            phi1(i,j) = phi(i,j)-dt*delta*sgn(phi0(i,j))*G(i,j);
            phi1(i,j)-= dt*(1-delta)*( sgn(phi0(i,j))*fabs(phi(i,j))-Dij )/1;
        }
    phi = phi1;
}

/*
 * check :
 *      retun |Grad(phi)|
 */

field LevelSet::check() const
{
    return (diff(_phi,0).square()+diff(_phi,1).square()).sqrt();
}

bool LevelSet::check(const double eps) const
{
    field mask = (5*1-_phi.abs()).max(0);
    mask /= mask+1e-15*1;

    std::cout << (mask*check()).abs().minCoeff() << " " << (mask*(check()-1).abs()).maxCoeff() << std::endl;

    if ( (mask*(check()-1).abs()).maxCoeff() > eps)
        return true;
    else
        return false;
}

/*
 * normal:
 *      return Grad(phi)/|Grad(phi)|
 */

std::vector<field> LevelSet::normal() const
{
    std::vector<field> n(2);
    for (int  i = 0; i<2; i++) n[i] = diff(_phi, i);
    field length = (n[0].square()+n[1].square()).sqrt();
    for (std::vector<field>::iterator it = n.begin(); it != n.end(); it++ ) *it/=length.max(1*1e-10);
    return n;
}

/*
 * curvature(vector<field3D> n):
 *      return div(n)
 */
field LevelSet::curvature(std::vector<field> &n) const
{
    field sigma = n[0]; sigma.setZero();
    for (int  i = 0; i<2; i++) {sigma += diff(n[i],i);}
    return sigma;
}

/*
 * curvature():
 *      return div(Grad(phi)/|Grad(phi)|)
 */
field LevelSet::curvature() const
{
    std::vector<field> n = this->normal();
    return curvature(n);
}

/*
 * Godunov()_v :
 *
 */
std::vector<std::vector<double>> LevelSet::godunov_v(const std::vector<std::vector<double>>& xm, const std::vector<std::vector<double>>& xp, const std::vector<std::vector<double>>& ym, const std::vector<std::vector<double>>& yp) const
{
    // field result = pos(_phi)*( (xm.max(0.).max(-xp.min(0.))).pow(2) + (ym.max(0.).max(-yp.min(0.))).pow(2)).sqrt();
    // result += neg(_phi)*( (xp.max(0.).max(-xm.min(0.))).pow(2) + (yp.max(0.).max(-ym.min(0.))).pow(2) ).sqrt();

    std::vector<std::vector<double>> result;
    result.resize(_phi.rows());
    for (int i=0;i<_phi.rows() ;i++) { result[i].resize(_phi.cols()); }

    double coef1( sqrt( pow(std::max( std::min(0.,-fmax(xp)), std::max(0., fmax(xm))), 2) + pow(std::max( std::min( 0., -fmin(yp)), std::max(0., fmax(ym))), 2)) );
    double coef2( sqrt( pow(std::max( std::min(0.,-fmax(xm)), std::max(0., fmax(xp))), 2) + pow(std::max( std::min( 0., -fmax(ym)), std::max(0., fmax(yp))), 2)) );

    for (int i; i<result.size(); i++)
    {
      for (int j; j<result[0].size(); j++)
      {
        result[i][j] = (pos_v(_phi_v,i,j) * coef1 + neg_v(_phi_v,i,j)*coef2) -1.;

      }
    }

    return result ; //- 1;
}

/*
 * Godunov() :
 *
 */
field LevelSet::godunov(const field& xm, const field& xp, const field& ym, const field& yp) const
{
    field result = pos(_phi)*( (xm.max(0.).max(-xp.min(0.))).pow(2) + (ym.max(0.).max(-yp.min(0.))).pow(2)).sqrt();
    result += neg(_phi)*( (xp.max(0.).max(-xm.min(0.))).pow(2) + (yp.max(0.).max(-ym.min(0.))).pow(2) ).sqrt();
    return result - 1;
}


/*
 * sgn(double X):
 *      return X/(|X|+eps)
 */
double LevelSet::sgn(double x) const
{
    return x/(fabs(x)+1*1e-15);
}

/*
 * sgn_sussman(field3D X):
 *      return X/sqrt(X^2+dx^2)
 */
field LevelSet::sgn_sussman(const field &X) const
{
    return X/(X.square()+pow(1,2)).sqrt();
}

/*
 * sgn_peng():
 *      return phi/sqrt(phi^2+|grad(phi)|^2+dx^2)
 */
field LevelSet::sgn_peng() const
{
    return _phi/(_phi.square()+this->check().square()*pow(1,2)).sqrt();
}

double LevelSet::boundary(const double ij, const double imj,  const double ipj, const double ijm, const double ijp) const
{
    double eps = 1*1e-15;

    double a = ( 1+std::min(ij*imj,0.0)/(fabs(ij*imj)+eps) )*( 1+std::min(ij*ipj,0.0)/(fabs(ij*ipj)+eps) );
    double b = ( 1+std::min(ij*ijm,0.0)/(fabs(ij*ijm)+eps) )*( 1+std::min(ij*ijp,0.0)/(fabs(ij*ijp)+eps) );

    return a*b;
}

/*
 * pos(X)_v :
 *      return max(X,0)/(|X|+eps)
 */
double LevelSet::pos_v(const std::vector<std::vector<double>>& X, int a, int b) const
{
    double max=std::max(fmax(X),0.);
    std::vector<std::vector<double>> newX;
    newX.resize(X.size());
    for (int i=0;i<newX.size() ;i++) { newX[i].resize(X[0].size()); }

    for (int i; i<newX.size(); i++)
    {
      for (int j; j<newX[0].size(); j++)
      {
        newX[i][j]= max/(abs(X[i][j])+1e-15);
      }
    }
    return newX[a][b];
}

/*
 * neg(X)_v :
 *      return -min(X,0)/(|X|+eps)
 */
double LevelSet::neg_v(const std::vector<std::vector<double>>& X, int a, int b) const
{
  double min=std::min(fmin(X),0.);
  std::vector<std::vector<double>> newX;
  newX.resize(X.size());
  for (int i=0;i<newX.size() ;i++) { newX[i].resize(X[0].size()); }

  for (int i; i<newX.size(); i++)
  {
    for (int j; j<newX[0].size(); j++)
    {
      newX[i][j]= min/(abs(X[i][j])+1e-15);
    }
  }
  return newX[a][b];
}

/*
 * pos(X) :
 *      return max(X,0)/(|X|+eps)
 */
field LevelSet::pos(const field& X) const
{
    return X.max(0.0)/(X.abs()+1*1e-15);
}

/*
 * neg(X) :
 *      return -min(X,0)/(|X|+eps)
 */
field LevelSet::neg(const field& X) const
{
    return -X.min(0.0)/(X.abs()+1*1e-15);
}

field LevelSet::diff(const field &F, bool index) const
{
  auto diffF = F; diffF.setZero();
  int N1(F.rows()), N2(F.cols()), d1(1.0), d2(1.0);
  if (!index) { // x
    diffF.block(1, 0, N1 - 2, N2) =
        0.5 / d1 * (F.block(2, 0, N1 - 2, N2) - F.block(0, 0, N1 - 2, N2));
    diffF.block(0, 0, 1, N2) =
        0.5 / d1 * (4.0 * F.block(1, 0, 1, N2) - F.block(2, 0, 1, N2) -
                    3.0 * F.block(0, 0, 1, N2));
    diffF.block(N1 - 1, 0, 1, N2) =
        0.5 / d1 *
        (-4.0 * F.block(N1 - 2, 0, 1, N2) + F.block(N1 - 3, 0, 1, N2) +
         3.0 * F.block(N1 - 1, 0, 1, N2));
  } else { // y
    diffF.block(0, 1, N1, N2 - 2) =
        0.5 / d2 * (F.block(0, 2, N1, N2 - 2) - F.block(0, 0, N1, N2 - 2));
    diffF.block(0, 0, N1, 1) =
        0.5 / d2 * (4.0 * F.block(0, 1, N1, 1) - F.block(0, 2, N1, 1) -
                    3.0 * F.block(0, 0, N1, 1));
    diffF.block(0, N2 - 1, N1, 1) =
        0.5 / d2 *
        (-4.0 * F.block(0, N2 - 2, N1, 1) + F.block(0, N2 - 3, N1, 1) +
         3.0 * F.block(0, N2 - 1, N1, 1));
  }
  return diffF;
}


std::vector<std::vector<double>> LevelSet::diff_v(const std::vector<std::vector<double>> &F, bool index) const
{
  auto diffF = F;

  for (int i=0; i<diffF.size(); i++)
  {
    for(int j=0; j<diffF[0].size(); j++ )
    {
      diffF[i][j] = 0.;
    }
  }
  int N1(F.size()), N2(F[0].size()), d1(1.0), d2(1.0);
  if (!index) { // x
    //diffF.block(1, 0, N1 - 2, N2) =
    //    0.5 / d1 * (F.block(2, 0, N1 - 2, N2) - F.block(0, 0, N1 - 2, N2));
    for (int i=1; i<N1-1; i++)
    {
      for (int j=0; j<N2; j++)
      {
        diffF[i][j] = 0.5 / d1 * (F[i+1][j] - F[i-1][j]);
      }
    }
    //diffF.block(0, 0, 1, N2) =
    //    0.5 / d1 * (4.0 * F.block(1, 0, 1, N2) - F.block(2, 0, 1, N2) -
    //                3.0 * F.block(0, 0, 1, N2));
    for (int j=0; j<N2; j++)
    {
      diffF[0][j] = 0.5 / d1 * (4.0 * F[1][j] - F[2][j] - 3.0 * F[0][j]);
    }
    //diffF.block(N1 - 1, 0, 1, N2) =
    //    0.5 / d1 *
    //    (-4.0 * F.block(N1 - 2, 0, 1, N2) + F.block(N1 - 3, 0, 1, N2) +
    //     3.0 * F.block(N1 - 1, 0, 1, N2));
    for (int j=0; j<N2; j++)
    {
      diffF[N1-1][j] = 0.5 / d1 * (-4.0 * F[N1-2][j] + F[N1-3][j] + 3.0 * F[N1-1][j]);
    }
  } else { // y
    //diffF.block(0, 1, N1, N2 - 2) =
    //    0.5 / d2 * (F.block(0, 2, N1, N2 - 2) - F.block(0, 0, N1, N2 - 2));
    for (int i=0; i<N1; i++)
    {
      for (int j=1; j<N2-1; j++)
      {
        diffF[i][j] = 0.5 / d2 * (F[i][j+1] - F[i][j-1]);
      }
    }
    //diffF.block(0, 0, N1, 1) =
    //    0.5 / d2 * (4.0 * F.block(0, 1, N1, 1) - F.block(0, 2, N1, 1) -
    //                3.0 * F.block(0, 0, N1, 1));
    for (int i=0; i<N1; i++)
    {
      diffF[i][0] = 0.5 / d2 * (4.0 * F[i][1] - F[i][2] - 3.0 * F[i][0]);
    }
    //diffF.block(0, N2 - 1, N1, 1) =
    //    0.5 / d2 *
    //    (-4.0 * F.block(0, N2 - 2, N1, 1) + F.block(0, N2 - 3, N1, 1) +
    //     3.0 * F.block(0, N2 - 1, N1, 1));
    for (int i=0; i<N1; i++)
    {
      diffF[i][N2-1] = 0.5 / d2 * (-4.0 * F[i][N2-2] + F[i][N2-3] + 3.0 * F[i][N2-1]);
    }
  }
  return diffF;
}

/*
 * diff_m_v:
 *      dérivé à gauche (m <-> moins)
 *
 */
std::vector<std::vector<double>> LevelSet::diff_m_v(const std::vector<std::vector<double>>& F, int index) const
{
    std::vector<std::vector<double>> diff_v = F;
    for (int i=0; i<diff_v.size(); i++)
    {
      for(int j=0; j<diff_v[0].size(); j++ )
      {
        diff_v[i][j] = 0.;
      }
    }
    for (int  i=0; i<F.size(); ++i)
        for (int  j=0; j<F[0].size(); ++j)
            diff_v[i][j] = diff_m_v(F, index, i, j);
    return diff_v;
}

double LevelSet::diff_m_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const {
    if(index==0) { //  x
        if(x>0)
            return 1.0/1 * (F[x][y]-F[x-1][y]);
        else
            return 0.0; // 1.0/1 * (F(x+1,y,z)-F(x,y,z));
    }
    else if (index==1) { // y
        if(y>0)
            return 1.0/1.0 * (F[x][y]-F[x][y-1]);
        else
            return 0.0; // 1.0/1.0 * (F(x,y+1,z)-F(x,y,z));
    }
    // NOT REACHED
    return 0.0;
}

/*
 * diff_p_v:
 *      dérivé à droite (p <-> plus)
 *
 */
std::vector<std::vector<double>> LevelSet::diff_p_v(const std::vector<std::vector<double>>& F, int index) const
{
    std::vector<std::vector<double>> diff_v = F;
    for (int i=0; i<diff_v.size(); i++)
    {
      for(int j=0; j<diff_v[0].size(); j++ )
      {
        diff_v[i][j] = 0.;
      }
    }

    for (int  i=0; i<F.size(); ++i)
        for (int  j=0; j<F[0].size(); ++j)
            diff_v[i][j] = diff_p_v(F, index, i, j);
    return diff_v;
}

double LevelSet::diff_p_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const
{
    if (index==0) { //  x
        if(x<F.size()-1)
            return 1.0/1 * (F[x+1][y]-F[x][y]);
        else
            return 0.0; // 1.0/1 * (F(x,y,z)-F(x-1,y,z));
    }
    else if (index==1){ // y
        if(y<F[0].size()-1)
            return 1.0/1.0 * (F[x][y+1]-F[x][y]);
        else
            return 0.0; // 1.0/1.0 * (F(x,y,z)-F(x,y-1,z));

    }
    // NOT REACHED
    return 0.0;
}

//fonction max min
// double fmax(std::vector<double> v)
// {
//   double max_v(v[0]);
//   for (size_t i = 1; i < v.size(); i++) {
//     if ( v[i]> max_v) {max_v = v[i];};
//   }
//   return max_v;
// }

double fmax(const std::vector<std::vector<double>> v)
{
  double max_v(v[0][0]);
  for (size_t i = 1; i < v.size(); i++) {
    for (size_t j = 0; j < v[0].size(); j++) {
      if ( v[i][j]> max_v) {max_v = v[i][j];};
    }
  }
  return max_v;
}

// double fmin(std::vector<double> v)
// {
//   double min_v(v[0]);
//   for (size_t i = 1; i < v.size(); i++) {
//     if ( v[i]< min_v) {min_v = v[i];};
//   }
//   return min_v;
// }

double fmin(std::vector<std::vector<double>> v)
{
  double min_v(v[0][0]);
  for (size_t i = 1; i < v.size(); i++) {
    for (size_t j = 0; j < v[0].size(); j++) {
      if ( v[i][j]< min_v) {min_v = v[i][j];};
    }
  }
  return min_v;
}

/*
 * diff_m:
 *      dérivé à gauche (m <-> moins)
 *
 */
field LevelSet::diff_m(const field& F, int index) const
{
    field diff = F; diff.setZero();
    for (int  i=0; i<F.rows(); ++i)
        for (int  j=0; j<F.cols(); ++j)
            diff(i,j) = diff_m(F, index, i, j);
    return diff;
}

double LevelSet::diff_m(const field& F, int index, int x, int y) const {
    if(index==0) { //  x
        if(x>0)
            return 1.0/1 * ( F(x,y)-F(x-1,y));
        else
            return 0.0; // 1.0/1 * (F(x+1,y,z)-F(x,y,z));
    }
    else if (index==1) { // y
        if(y>0)
            return 1.0/1.0 * (F(x,y)-F(x,y-1));
        else
            return 0.0; // 1.0/1.0 * (F(x,y+1,z)-F(x,y,z));
    }
    // NOT REACHED
    return 0.0;
}

/*
 * diff_p:
 *      dérivé à droite (p <-> plus)
 *
 */
field LevelSet::diff_p(const field& F, int index) const
{
    field diff = F; diff.setZero();
    for (int  i=0; i<F.rows(); ++i)
        for (int  j=0; j<F.cols(); ++j)
            diff(i,j) = diff_p(F, index, i, j);
    return diff;
}

double LevelSet::diff_p(const field& F, int index, int x, int y) const
{
    if (index==0) { //  x
        if(x<F.rows()-1)
            return 1.0/1 * ( F(x+1,y)-F(x,y));
        else
            return 0.0; // 1.0/1 * (F(x,y,z)-F(x-1,y,z));
    }
    else if (index==1){ // y
        if(y<F.cols()-1)
            return 1.0/1.0 * (F(x,y+1)-F(x,y));
        else
            return 0.0; // 1.0/1.0 * (F(x,y,z)-F(x,y-1,z));

    }
    // NOT REACHED
    return 0.0;
}



#endif
