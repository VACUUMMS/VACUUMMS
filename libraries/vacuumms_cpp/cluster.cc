/* vacuumms/cluster.cc */

#include <vacuumms/cavity.hh>


void CavityCluster::pushBack(Cavity cavity)
{   
    configuration.pushBack(cavity);
}


std::vector<CavityCluster> CavityConfiguration::generateClusters()
{
    std::vector<CavityCluster> clusters;
    
    // check non-zero box dimensions
    
    if (
        (box_dimensions[0] <= 0.0) || 
        (box_dimensions[1] <= 0.0) || 
        (box_dimensions[2] <= 0.0) 
        )
    {
        std::cout << "One or more box dimesions is <= 0. Failing." << std::endl;
        std::cout << box_dimensions[0] << std::endl;
        std::cout << box_dimensions[1] << std::endl;
        std::cout << box_dimensions[2] << std::endl;
        return clusters;
    }

    // find all pairs
    
    std::vector<std::pair<int, int>> pairs;

    for (int i = 0; i < records.size() - 1; i++)
    for (int j = i + 1; j < records.size(); j++)
    {
        for (vacuumms_float shift_x = -box_dimensions[0]; shift_x <= box_dimensions[0]; shift_x += box_dimensions[0])
        for (vacuumms_float shift_y = -box_dimensions[1]; shift_y <= box_dimensions[1]; shift_y += box_dimensions[1])
        for (vacuumms_float shift_z = -box_dimensions[2]; shift_z <= box_dimensions[2]; shift_z += box_dimensions[2])
        {
            vacuumms_float dx2 = (shift_x + records[i].x - records[j].x) * (shift_x + records[i].x - records[j].x);
            vacuumms_float dy2 = (shift_y + records[i].y - records[j].y) * (shift_y + records[i].y - records[j].y);
            vacuumms_float dz2 = (shift_z + records[i].z - records[j].z) * (shift_z + records[i].z - records[j].z);

            if ((dx2 + dy2 + dz2) < ((records[i].d + records[j].d) * (records[i].d + records[j].d) * .25))
                pairs.emplace_back(i,j);
        }
    }

    // initialize clusters to index of cavities
    
    std::vector<int> cluster_numbers;
    for (int i = 0; i < records.size(); i++) cluster_numbers.emplace_back(i);

    // build clusters

    for (int i=0; i<records.size(); i++) 
    for (int j=0; j<pairs.size(); j++) 
        if (pairs[j].first == i) cluster_numbers[pairs[j].second] = cluster_numbers[i];

    // build cluster objects 
    // cluster_numbers are the one-to-one mapping of cavity number to cluster number
    
    for (int i = 0; i < cluster_numbers.size(); i++) // for each cluster number
    {
        CavityCluster cluster; // create an empty cluster

        // find all cavities mapped to that cluster number
        for (int j = 0; j < cluster_numbers.size(); j++) 
            if (cluster_numbers[j] == i) cluster.pushBack(records[j]);
        
        // size zero clusters get skipped and discarded
        if (cluster.getSize() > 0) clusters.push_back(cluster);
    }

    return clusters;

} // end generateClusters()


int CavityCluster::getSize()
{
    return configuration.getSize();
}


#ifdef BUILD_PYBIND_BINDINGS

pybind11::str CavityCluster::__repr__()
{
    return configuration.__repr__();
}

#endif

