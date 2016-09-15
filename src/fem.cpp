#include "fem.h"
#include <Eigen/SparseCholesky>
#include <iostream>

Triangle::Triangle(Vector a, Vector b, Vector c, std::array< int, int(3) > indices)
{
    this->points[0] = a;
    this->points[1] = b;
    this->points[2] = c;
    this->indices = indices;
    Eigen::Matrix<double, 2, 3> B_eta;
    Eigen::Matrix<double, 2, 2> J;
    this->center.setZero();
    for (auto point: points)
        this->center += point;
    this->center /= 3;
    
    B_eta << -1, 1, 0, -1, 0, 1;
    pos_mat << points[0].x(), points[0].y(),
               points[1].x(), points[1].y(),
               points[2].x(), points[2].y();    
    J = B_eta * pos_mat;
    this->area = J.determinant() / 2;
    this->B = J.inverse() * B_eta;
}

void Triangle::add_fem_equation(std::vector< Eigen::Triplet< double > >& K_g, Eigen::VectorXd & rhs_g)
{
    Eigen::Matrix<double, 3, 3> K_l;
    Eigen::Matrix<double, 3, 1> rhs_l;
    K_l = (this->B.transpose() * this->B) * this->area;
    rhs_l = this->area * this->B.transpose() * Vector(this->center.y(), - this->center.x());
    int row_l, col_l = 0;
    for (int row_g: this->indices)
    {
        if (row_g != 0)
        {
            col_l = 0;
            for (int col_g: this->indices)
            {
                K_g.push_back(Eigen::Triplet<double> (row_g, col_g, K_l(row_l, col_l)));
                col_l++;
            }
            rhs_g(row_g) += rhs_l(row_l);
        }
        row_l++;
    }
}

Vector Triangle::get_grad(const Eigen::VectorXd& sol)
{
    return B * Eigen::Matrix<double, 3, 1>(sol[indices[0]], sol[indices[1]], sol[indices[2]]);
}



TorsionFemCase::TorsionFemCase(std::vector<Vector> vertices, std::vector<std::array<int, 3>> triangles)
{
    for (auto tri_index: triangles)
    {
        this->triangles.push_back(
            Triangle(
                vertices[tri_index[0]],
                vertices[tri_index[1]],
                vertices[tri_index[2]],
                tri_index));
    }
    
    this->mat_size = vertices.size();
}

void TorsionFemCase::run()
{
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    Eigen::VectorXd rhs(mat_size);
    rhs.setZero();
    std::vector<Eigen::Triplet<double>> mat_entries;
    for (auto triangle: this->triangles)
    {
        triangle.add_fem_equation(mat_entries, rhs);
    }
//     fixing w at one position
    mat_entries.push_back(Eigen::Triplet<double>(1, 0, 0));
    rhs[0] = 0;
    mat.setFromTriplets(mat_entries.begin(), mat_entries.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(mat);
    this->sol = cg.solve(rhs);
}

double TorsionFemCase::get_polar_moment()
{
    double polar_moment = 0;
    for (auto triangle : this->triangles)
        polar_moment += triangle.center.dot(triangle.center) * triangle.area;
    return polar_moment;
}

double TorsionFemCase::get_torsion_moment()
{
    double torsion_moment = 0;
    for (auto triangle : this->triangles)
    {
        Vector grad_w = triangle.get_grad(this->sol);
        torsion_moment += triangle.center.dot(triangle.center) * triangle.area;
        torsion_moment -= triangle.center.y() * grad_w.x() * triangle.area;
        torsion_moment += triangle.center.x() * grad_w.y() * triangle.area;
    }
    return torsion_moment;
}

std::vector< std::array< double, 2> > TorsionFemCase::get_stress()
{
    std::vector< std::array< double, 2 > > stress;
    for (auto triangle : this->triangles)
    {
        Vector grad_w = triangle.get_grad(this->sol);

        std::array<double, 2> temp;
        temp[0] = - triangle.center.y() + grad_w.x();
        temp[1] =   triangle.center.x() + grad_w.y();
        stress.push_back(temp);
    }
    return stress;
}

std::vector< double > TorsionFemCase::get_w()
{
    std::vector< double > w;
    for (int i=0; i<this->sol.size(); i++)
    {
        w.push_back(this->sol(i));
    }
    return w;
}

