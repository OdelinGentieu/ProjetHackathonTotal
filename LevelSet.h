#ifndef LEVELSET_H
#define LEVELSET_H

#include "Dense"
#include "Image.h"
#include <iostream>

#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class LevelSet
{
private:
    field& _phi;
    std::vector<std::vector<double>>& _phi_v;

    std::vector<field> normal() const;
    field curvature(std::vector<field> &n) const;

    double sgn(double x) const;
    field sgn_peng() const;
    field sgn_sussman(const field &X) const;

    double pos_v(const std::vector<std::vector<double>>& X, int a, int b) const;
    field pos(const field& X) const;
    double neg_v(const std::vector<std::vector<double>>& X, int a, int b) const;
    field neg(const field& X) const;

    field godunov(const field& xm, const field& xp, const field& ym, const field& yp) const;
    std::vector<std::vector<double>> godunov_v(const std::vector<std::vector<double>>& xm, const std::vector<std::vector<double>>& xp, const std::vector<std::vector<double>>& ym, const std::vector<std::vector<double>>& yp) const;

    double boundary(const double ij, const double imj,  const double ipj, const double ijm, const double ijp) const;

    void firstOrderUpwind_v(std::vector<std::vector<double>>& phi, const std::vector<std::vector<double>>& phi0, const double dt);
    void firstOrderUpwind(field& phi, const field& phi0, const double dt);


    double diff_m_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const;
    double diff_m(const field& F, int index, int x, int y) const;

    double diff_p_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const;
    double diff_p(const field& F, int index, int x, int y) const;

    std::vector<std::vector<double>> diff_m_v(const std::vector<std::vector<double>>& F, int index) const;
    field diff_m(const field& F, int index) const;

    std::vector<std::vector<double>> diff_p_v(const std::vector<std::vector<double>>& F, int index) const;
    field diff_p(const field& F, int index) const;

    std::vector<std::vector<double>> diff_v(const std::vector<std::vector<double>>&F, bool index) const;
    field diff(const field &F, bool index) const;


  public:
      LevelSet(field& phi);//:_phi(phi){}
      ~LevelSet() = default;
      void redistancing(const int max_iter);
      void redistancing(const double eps);

      field curvature() const;
      field check() const;
      bool check(const double eps) const;
      
      double fmax(const std::vector<std::vector<double>> v) const;
      double fmin(const std::vector<std::vector<double>> v) const;
};

#endif // LEVELSET_H
