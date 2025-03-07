
#include <iostream>
#include <vector>
#include <voro++.hh>

#include <vacuumms/voronoi.hh>

//using namespace voro;


Voronoi::Voronoi()
{}

Voronoi::Voronoi(Configuration gfg, Parameters prm)
{
    // Step 1: Define a container (non-periodic cube from -1 to 1)
    // con = container(-1, 1, -1, 1, -1, 1, 6, 6, 6, false, false, false, 8);
    int block_param = 6;
    container = std::make_unique<voro::container>(0, gfg.box_x, 0, gfg.box_y, 0, gfg.box_z, block_param, block_param, block_param, true, true, true, 8);

/*    // Step 2: Add particles
    con.put(0, 0.0, 0.0, 0.0); // Particle at origin
    con.put(1, 0.5, 0.5, 0.5); // Another particle
*/
    for (int i = 0; i<gfg.getSize(); i++)
    {
        container->put(i, gfg.records[i].x, gfg.records[i].y, gfg.records[i].z);
    }

    // Step 3: Loop over all cells and compute geometry
    voro::voronoicell c;
    voro::c_loop_all cl(*container);
    if (cl.start()) do 
    {
        if (container->compute_cell(c, cl)) {
            std::cout << "Particle " << cl.pid() << ":\n";

            // Get vertices
            std::vector<double> vertices;
            c.vertices(vertices); // Fills with x, y, z for each vertex

            std::cout << "Vertices:\n";
            for (size_t i = 0; i < vertices.size(); i += 3) {
                std::cout << i/3 << ": (" << vertices[i] << ", " << vertices[i+1] 
                          << ", " << vertices[i+2] << ")\n";
            }

            // Get face vertex indices
            std::vector<int> face_vertices;
            c.face_vertices(face_vertices); // Fills with face data: [n, v1, v2, ..., vn, n, ...]

            // Extract edges from face vertices
            std::vector<VoronoiEdge> edges;
            int pos = 0;
            while (pos < face_vertices.size()) {
                int n = face_vertices[pos]; // Number of vertices in this face
                for (int i = 0; i < n; i++) {
                    int v1 = face_vertices[pos + 1 + i];
                    int v2 = face_vertices[pos + 1 + (i + 1) % n]; // Wrap around to first vertex
                    VoronoiEdge e(v1, v2);
                    // Avoid duplicates (since edges are shared between faces)
                    if (std::find(edges.begin(), edges.end(), e) == edges.end()) {
                        edges.push_back(e);
                    }
                }
                pos += n + 1; // Move to next face
            }

            std::cout << "VoronoiEdges:\n";
            for (const VoronoiEdge& e : edges) {
                std::cout << "(" << e.v1 << ", " << e.v2 << ")\n";
            }
        }
    } while (cl.inc());

}

void Voronoi::execute()
{}
