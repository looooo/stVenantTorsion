#include "bem.h"
#include <vector>
#include <tuple>
#include <memory>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include <iostream>


DirichletPanel::DirichletPanel(VectorPtr p1, VectorPtr p2): Panel(p1, p2)
{
    bc_neumann = - this->normal.x() * this->center.y() +
                   this->normal.y() * this->center.x();
}

std::tuple<double, double> DirichletPanel::get_influence(Vector target)
{
    return std::make_tuple(
        dipole_0(target, *this), 
        monopole_0(target, *this)* this->bc_neumann);
}

std::tuple<Vector, Vector> DirichletPanel::get_grad_influence(Vector target)
{
    return std::make_tuple(
        dipole_0_v(target, *this), 
        monopole_0_v(target, *this)* this->bc_neumann);
}

void DirichletPanel::add_boundary_equation(
    std::vector<DPanelPtr> panels, Eigen::MatrixXd& mat, Eigen::VectorXd& rhs, int i)
{
    std::tuple<double, double> infl;
    for (int j=0; j < panels.size(); j++)
    {
            infl = panels[j]->get_influence(this->center);
            mat(i, j) = std::get<0>(infl);
            rhs(i) += std::get<1>(infl);
    }
}

double DirichletPanel::get_polar_influence()
{
    return (this->center.x() * pow(this->center.y(), 2) * this->normal.x() +
            this->center.y() * pow(this->center.x(), 2) * this->normal.y()) * this->area;
}

TorsionBemCase::TorsionBemCase(std::vector< DPanelPtr > panels)
{
    this->panels = panels;
}

void TorsionBemCase::run()
{
    Eigen::MatrixXd mat(this->panels.size(), this->panels.size());
    Eigen::VectorXd rhs(this->panels.size());
    this->sol.resize(this->panels.size());
    mat.setZero();
    rhs.setZero();
    int i = 0;
    for(auto panel: this->panels)
        panel->add_boundary_equation(this->panels, mat, rhs, i++);
    // fixing w at one position
    mat.row(0).setZero();
    mat(0, 0) = 1;
    rhs(0) = 0;
    this->sol = mat.lu().solve(rhs);
}

double TorsionBemCase::get_w(Vector point)
{
    std::tuple<double, double> infl;
    double t = 0;
    int i = 0;
    
    for (auto panel: this->panels)
    {
        infl = panel->get_influence(point);
        t += std::get<0>(infl) * this->sol[i] - std::get<1>(infl);
        i++;
    }
    return t;
}

Vector TorsionBemCase::get_stress(Vector point)
{
    std::tuple<Vector, Vector> infl;
    Vector q(0, 0);
    int i = 0;
    
    for (auto panel: this->panels)
    {
        infl = panel->get_grad_influence(point);
        q += std::get<0>(infl) * this->sol[i] - std::get<1>(infl);
        return Vector(q.x() - point.y(), q.y() + point.x());
        i++;
    }
}

double TorsionBemCase::get_polar_moment()
{
    double pol = 0;
    for (auto panel: this->panels)
        pol += panel->get_polar_influence();
    return pol;
}

double TorsionBemCase::get_torsion_moment()
{
    double s = 0;
    int i = 0;
    for (auto panel: this->panels)
    {
        s += this->sol[i] * panel->bc_neumann * panel->area;
        i++;
    }
    return s + this->get_polar_moment();
}