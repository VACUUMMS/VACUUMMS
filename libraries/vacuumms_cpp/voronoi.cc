/* libraries/vacuumms_cpp/voronoi.cc */

#include <iostream>
#include <vector>
#include <voro++.hh>

#include <vacuumms/voronoi.hh>

//using namespace voro;


Voronoi::Voronoi()
{}

Voronoi::Voronoi(Configuration gfg, Parameters prm)
{
    int block_param = 6;
	// class member container
    container = std::make_unique<voro::container>(0, gfg.box_x, 0, gfg.box_y, 0, gfg.box_z, block_param, block_param, block_param, true, true, true, 8);

    // Step 2: Add particles (example coordinates from your POV-Ray spheres, simplified)
    for (int i = 0; i<gfg.getSize(); i++)
    {
        container->put(i, gfg.records[i].x, gfg.records[i].y, gfg.records[i].z);
    }

    // Step 3: Extract global vertices and edges
    //    std::vector<VoronoiVertex> global_vertices;
    std::map<VoronoiVertex, int> vertex_index_map; // Map vertices to global indices
    std::set<VoronoiEdge> global_edges;            // Set for unique edges

    voro::voronoicell c;
    voro::c_loop_all cl(*container);

    if (cl.start()) {
        do {
            if (container->compute_cell(c, cl)) {
                // Local vertices for this cell
                std::vector<double> vertex_coords;
                c.vertices(vertex_coords);

                // Map local vertex indices to global indices
                std::vector<int> local_to_global(vertex_coords.size() / 3);
                for (size_t i = 0; i < vertex_coords.size(); i += 3) {
                    int idx = i/3;
                    VoronoiVertex v(i, vertex_coords[i], vertex_coords[i + 1], vertex_coords[i + 2]);
                    auto it = vertex_index_map.find(v);
                    if (it == vertex_index_map.end()) {
                        int global_idx = global_vertices.size();
                        global_vertices.push_back(v);
                        vertex_index_map[v] = global_idx;
                        local_to_global[i / 3] = global_idx;
                    } else {
                        local_to_global[i / 3] = it->second;
                    }
                }

                // Extract face vertex indices and compute edges
                std::vector<int> face_vertices;
                c.face_vertices(face_vertices);

                int pos = 0;
                while (pos < face_vertices.size()) {
                    int n = face_vertices[pos]; // Number of vertices in this face
                    for (int i = 0; i < n; ++i) {
                        int local_v1 = face_vertices[pos + 1 + i];
                        int local_v2 = face_vertices[pos + 1 + (i + 1) % n]; // Wrap around
                        int global_v1 = local_to_global[local_v1];
                        int global_v2 = local_to_global[local_v2];
                        global_edges.emplace(global_v1, global_v2);
                    }
                    pos += n + 1;
                }
            }
        } while (cl.inc());
    }

    // Convert global_edges set to vector for consistency
    //std::vector<VoronoiEdge> edge_list(global_edges.begin(), global_edges.end());
    edge_list = std::vector<VoronoiEdge> (global_edges.begin(), global_edges.end());

    // Step 4: Output the global data
//    std::cout << "\nGlobal VoronoiVertices (" << global_vertices.size() << "):\n";
//    std::cout << global_vertices.size();

    std::cout << "\nGlobal Unique VoronoiEdges (" << edge_list.size() << "):\n";
	std::cout << edge_list.size();

}


void Voronoi::execute()
{}

std::vector<VoronoiVertex> Voronoi::getVertices()
{
    return global_vertices;
}

std::vector<VoronoiEdge> Voronoi::getEdges()
{
    return edge_list;
}

#ifdef BUILD_PYBIND_BINDINGS
	
pybind11::str VoronoiVertex::__repr__()
{
    pybind11::str retval("(");
    retval = retval +
        pybind11::str(std::to_string(index)) +
        pybind11::str(": ") +
        pybind11::str(std::to_string(x)) +
        pybind11::str(", ") +
        pybind11::str(std::to_string(y)) +
        pybind11::str(", ") +
        pybind11::str(std::to_string(z)) +
        pybind11::str(")\n");
    return retval;
}

	
pybind11::str VoronoiEdge::__repr__()
{
    pybind11::str retval("(");
    retval = retval +
        pybind11::str(std::to_string(v1)) +
        pybind11::str(", ") +
        pybind11::str(std::to_string(v2)) +
        pybind11::str(")\n");
    return retval;
}

#endif

