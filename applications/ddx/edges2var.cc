/* edges2pairs.cc */

//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

#include <vacuumms/cavity.hh>
#include <vacuumms/vertex.hh>
#include <vacuumms/edge.hh>
#include <vacuumms/pair.hh>

#include <vacuumms/param.hh>
#include <vacuumms/types.h>

double box_x=0.0, box_y=0.0, box_z=0.0;
char *verts_filename;
char *edges_filename;
char *map_filename;

int main(int argc, char *argv[])
{
    // Parse command line options

    setCommandLineParameters(argc, argv);
    if (getFlagParam((char*)"-usage"))
    {
        printf("\n");
        printf("vacuumms/edges2var: \n");
        printf("\n");
        printf("edges2var takes an indexed list of unique cavities (vcv format),\n");
        printf("a list of edges (vertex pairs) from the voronoi diagram, and a\n");
        printf("map file that maps vertices by index to the cavities generated by\n");
        printf("insertion at each vertex. A vertices file is also needed to map\n");
        printf("the vertices to the edges. This generates a list of **pairs** of\n");
        printf("cavities based on the voronoi initial condition that can be used\n");
        printf("as inputs to a variational iteration of these pair paths, giving\n");
        printf("an image of likely paths through the material.\n");
        printf("\n");
        printf("usage:       edges2var");
        printf("             -box [ n.nn n.nn n.nn ] \n");
        printf("             -edges_file [ edges.edg ] \n");
        printf("             -verts_file [ vertices.vrt ] \n");
        printf("             -map_file [ mapfile.map ] \n");
        printf("             < cavities.vnq (read formatted input from stdin) \n");
        printf("             [ where .vnq record format is:\n");
        printf("             %%d     %%f     %%f     %%f     %%f     %%f\n");
        printf("             index     x     y     z     d     drift ]\n"); 
        printf("\n");
        return 0;
    }

    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    if (box_x * box_y * box_z <= 0.0) 
    {
        fprintf(stderr, "vacuumms/edges2pairs: please specify box size.\n");
        exit(2);
    }

    getStringParam((char *)"-verts_file", &verts_filename);
    getStringParam((char *)"-map_file", &map_filename);
    getStringParam((char *)"-edges_file", &edges_filename);

    if (verts_filename == NULL || edges_filename == NULL || map_filename == NULL) 
    {
        fprintf(stderr, "vacuumms/edges2pairs: please specify verts_file, edges_file and map_file.\n");
        exit(2);
    }

    VertexList verts(verts_filename);

    EdgeList edges(edges_filename, verts);

    // we've got vertices and edges lists. Now we refine:

    // load the vcav output, so we can map verts to cavities. 

    CavityConfiguration vnq;

    int index;
    vacuumms_float x, y, z, d, drift;

    while (!feof(stdin))
    {
        fscanf(stdin, "%d%f\t%f\t%f\t%f%f\n", &index, &x, &y, &z, &d, &drift);
        int record_number = vnq.getSize();
        Cavity cavity(index, x, y, z, d, drift);
        // cavity.setForeignKey(??);
        vnq.records.push_back(cavity);
        // vnq.recordAt(record_number).setForeignKey(index);
    }

    // vertex-to-cavity map
    IndexPairList v2c_map(map_filename);

    // everything loaded in, so now start culling pairs for output

    IndexPairList output_list;

    // Loop over the edges. For each edge, turn the vertex index pair into a cavity index pair
    for (int i = 0; i < edges.getSize(); i++)
    {
        Edge edge = edges.recordAt(i);

        int A_found = 0;
        int B_found = 0;

        int A_out = -1;
        int B_out = -1;

        for (int j = 0; j < v2c_map.getSize(); j++)
        { 
            IndexPair v2c_ip = v2c_map.recordAt(j);
            
            // A
            if (!A_found && edge.A == v2c_ip.A) 
            {
                A_found = 1;
                A_out = v2c_ip.B;
            }

            // B
            if (!B_found && edge.B == v2c_ip.A) 
            {
                B_found = 1;
                B_out = v2c_ip.B;
            }

            if (A_found && B_found) 
            {
                break;
            }
        }
    
        if (!A_found || !B_found)
        {
            fprintf(stderr, "vacuumms/edges2pairs: failed to find A(%d) or B(%d)\n", A_found, B_found);
            exit(2);
        }

        // push the record onto the list
        output_list.pushBack(IndexPair(A_out, B_out));
    }

    // once the list is complete, remove duplicates and self-pairs

    for (int i=0;i<output_list.getSize(); i++) 
    {
        IndexPair ip_i = output_list.recordAt(i);
        for (int j=i+1; j<output_list.getSize(); j++) 
        {
            IndexPair ip_j = output_list.recordAt(j);
            // check forward and back
            if 
            (
                ((ip_i.A == ip_j.A) && (ip_i.B == ip_j.B)) || 
                ((ip_i.A == ip_j.B) && (ip_i.B == ip_j.A)) 
            ) 
            {
                // duplicate found 
                output_list.deleteRecordAt(j--); // adjust index j for deleted record
            }
        }
    }

    for (int i=0;i<output_list.getSize(); i++) 
    {
        IndexPair ip_i = output_list.recordAt(i);
        if (ip_i.A == ip_i.B) // check self
        {
            output_list.deleteRecordAt(i--); // adjust index i for deleted record
        }
    }

    // then turn this list/map into start/end points based on the cavs indexed

    for (int i=0;i<output_list.getSize(); i++) 
    {
        IndexPair ip_i = output_list.recordAt(i);
        vacuumms_float xA, yA, zA, xB, yB, zB;

        int found_A = 0;
        int found_B = 0;
        
        // map the A part of the output record
        for (int j=0; j < vnq.getSize(); j++)
        {
            if (ip_i.A == vnq.recordAt(j).index)
            {
                xA = vnq.recordAt(j).x;
                yA = vnq.recordAt(j).y;
                zA = vnq.recordAt(j).z;
                found_A = 1;
                break;
            }
        }

        // map the B part of the output record
        for (int j=0; j < vnq.getSize(); j++)
        {
            if (ip_i.B == vnq.recordAt(j).index)
            {
                xB = vnq.recordAt(j).x;
                yB = vnq.recordAt(j).y;
                zB = vnq.recordAt(j).z;
                found_B = 1;
                break;
            }
        }

        // Generate the output record

        printf("%f\t%f\t%f\t%f\t%f\t%f\n", xA, yA, zA, xB, yB, zB);
    }

    return 0;
}

