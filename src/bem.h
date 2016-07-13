#ifndef BEM
#define BEM


#include "laplace/laplaceKern2D.h"

using namespace laplaceKern2D;


class DirichletPanel;
typedef DirichletPanel DPanel;
typedef std::shared_ptr<DPanel> DPanelPtr;

struct DirichletPanel : Panel
{
    DirichletPanel(VectorPtr, VectorPtr);
    double bc_neumann;
    std::tuple<double, double> get_influence(Vector);
    std::tuple<Vector, Vector> get_grad_influence(Vector);
    void add_boundary_equation(std::vector<DPanelPtr>, Eigen::MatrixXd&, Eigen::VectorXd&, int);
    double get_polar_influence();
};

struct TorsionBemCase
{
    TorsionBemCase(std::vector<DPanelPtr>);
    std::vector<DPanelPtr> panels;
    Eigen::VectorXd sol;
    void run();
    double get_polar_moment();
    double get_torsion_moment();
    double get_w(std::array<double, 2> point);
    Vector get_stress(std::array<double, 2> point);
};

#endif