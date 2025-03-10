/* vacuumms/voronoi.hh */

#pragma once

#include <voro++.hh>

#include <vacuumms/exports.hh>
#include <vacuumms/types.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>

class VoronoiVertex
{
    public: 

    VoronoiVertex(){}
    VoronoiVertex(int idx, vacuumms_float x_, vacuumms_float y_, vacuumms_float z_) : index(idx), x(x_), y(y_), z(z_) {}
    int index;
    vacuumms_float x;
    vacuumms_float y;
    vacuumms_float z;

    bool operator<(const VoronoiVertex& other) const {
        return x < other.x || (x == other.x && (y < other.y || (y == other.y && z < other.z)));
    }

#ifdef BUILD_PYBIND_BINDINGS
    pybind11::str __repr__();    
#endif
    
};


class VoronoiEdge
{
    public: 

    VoronoiEdge(){}
    int v1, v2;
    VoronoiEdge(int v1_, int v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}
    bool operator<(const VoronoiEdge& other) const {
        return v1 < other.v1 || (v1 == other.v1 && v2 < other.v2);
    }

#ifdef BUILD_PYBIND_BINDINGS
    pybind11::str __repr__();    
#endif
    
};


class Voronoi : public Operation
{
    public:

        Voronoi();
        Voronoi(Configuration gfg, Parameters prm);
        std::vector<VoronoiVertex> getVertices();
        std::vector<VoronoiEdge> getEdges();
        void execute();

    private:
        
        Configuration cfg;
        Parameters prm;
        std::vector<VoronoiVertex> global_vertices;
        std::vector<VoronoiEdge> edge_list;
        std::unique_ptr<voro::container> container;
};
