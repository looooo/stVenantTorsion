#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Eigen/Core"
#include "bem.h"
#include "fem.h"

#include <vector>
#include <tuple>
#include <memory>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void bem_from_python(TorsionBemCase& cls, std::vector<std::vector<std::array<double, 2>>> vertices)
{
    
    std::vector<DPanelPtr> panels;
    for (auto boundary: vertices)
    {
        int vertex_index = 0;
        for (auto vertex: boundary)
        {
            if (vertex != boundary.back())
            {
                VectorPtr v1 (new Vector(vertex[0], vertex[1]));
                VectorPtr v2 (new Vector(boundary[vertex_index + 1][0], boundary[vertex_index + 1][1]));
                DPanelPtr panel (new DirichletPanel(v1, v2));
                panels.push_back(panel);
            }
            else
            {
                VectorPtr v1 (new Vector(vertex[0], vertex[1]));
                VectorPtr v2 (new Vector(boundary.front()[0], boundary.front()[1]));
                DPanelPtr panel (new DirichletPanel(v1, v2));
                panels.push_back(panel);
            }
            vertex_index ++;
        }
    }
    new (&cls) TorsionBemCase(panels);
}

void fem_from_python(TorsionFemCase& cls, 
                     std::vector<std::array<double, 2>> vertices,
                     std::vector<std::array<int, 3>> triangles)
{
    std::vector<Vector> vertices_new;
    for (auto vertex : vertices)
        vertices_new.push_back(Vector(vertex[0], vertex[1]));
    new (&cls) TorsionFemCase(vertices_new, triangles);
}

void init_stVenant(py::module &m)
{

    py::class_<TorsionBemCase>(m, "TorsionBemCase", "using boundary elements")
        .def("__init__", &bem_from_python)
        .def("get_polar_moment", &TorsionBemCase::get_polar_moment)
        .def("get_torsion_moment", &TorsionBemCase::get_torsion_moment)
        .def("get_w", &TorsionBemCase::get_w)
        .def("get_stress", &TorsionBemCase::get_stress)
        .def("run", &TorsionBemCase::run);

    py::class_<TorsionFemCase>(m, "TorsionFemCase", "using finite elements")
        .def("__init__", &fem_from_python)
        .def("get_polar_moment", &TorsionFemCase::get_polar_moment)
        .def("get_torsion_moment", &TorsionFemCase::get_torsion_moment)
        .def("get_w", &TorsionFemCase::get_w)
        .def("get_stress", &TorsionFemCase::get_stress)
        .def("run", &TorsionFemCase::run)
        .def_static("bc", &TorsionFemCase::bc);
};


PYBIND11_PLUGIN(stVenant)
{
    py::module m("stVenant", "torsion fem/bem");
    init_stVenant(m);
    return m.ptr();
}

