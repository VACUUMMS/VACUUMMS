/* vacuumms/voronoi.hh */

#include <voro++.hh>

#include <vacuumms/types.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>

struct VoronoiEdge
{
    int v1;
    int v2;

    VoronoiEdge(int a, int b) : v1(a), v2(b) {}

    // Edges are equal regardless of orientation
    bool operator==(const VoronoiEdge& other) const {
        return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
    }
};

struct VoronoiVertex
{
    int index;
    vacuumms_float x;
    vacuumms_float y;
    vacuumms_float z;
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
        std::vector<VoronoiVertex> Vertices;
        std::vector<VoronoiEdge> Edges;
        std::unique_ptr<voro::container> container;
};
