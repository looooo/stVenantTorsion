#ifndef FEM
#define FEM

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include <set>

typedef Eigen::Vector2d Vector;

class Triangle;

struct Triangle
{
    std::array<Vector,3> points;
    std::array<int, 3> indices;
    double area;
    Vector center;
    Eigen::Matrix<double, 2, 3> B;
    Eigen::Matrix<double, 3, 2> pos_mat;

    Triangle(Vector, Vector, Vector, std::array<int, 3>);
    Vector get_grad(const Eigen::VectorXd & sol);
    void add_fem_equation(std::vector<Eigen::Triplet< double>>& K_g, Eigen::VectorXd & rhs_g);
};


struct TorsionFemCase
{
    std::vector<Triangle> triangles;
    int mat_size;
    Eigen::VectorXd sol;
    TorsionFemCase(std::vector<Vector> vertices, std::vector<std::array<int, 3>> triangles);
    void run();
    double get_polar_moment();
    double get_torsion_moment();
    std::vector<double> get_w();
    std::vector<std::array<double, 2>> get_stress();
};



#endif