/* vacuumms/cluster.cc */

#include <vacuumms/cavity.hh>

/*
CavityCluster::CavityCluster()
{
}
*/


void CavityCluster::pushBack(Cavity cavity)
{   
    configuration.pushBack(cavity);
}


std::vector<CavityCluster> CavityConfiguration::generateClusters()
{
    std::vector<CavityCluster> clusters;
// need to replace these
// int number_of_clusters=0;
// int number_of_cavities=getSize();
// int number_of_pairs=0;

//double x[MAX_CAVITIES], y[MAX_CAVITIES], z[MAX_CAVITIES], d[MAX_CAVITIES];
//int cluster_number[MAX_CAVITIES];
//int cavityA[MAX_PAIRS], cavityB[MAX_PAIRS];

    // findAllPairs()
    
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

//            vacuumms_float dy2 = (shift_y + y[i] - y[j]) * (shift_y + y[i] - y[j]);
//            vacuumms_float dz2 = (shift_z + z[i] - z[j]) * (shift_z + z[i] - z[j]);
  
//            if ((dx2 + dy2 + dz2) < ((d[i] + d[j]) * (d[i] + d[j]) * .25))
            if ((dx2 + dy2 + dz2) < ((records[i].d + records[j].d) * (records[i].d + records[j].d) * .25))
                pairs.emplace_back(i,j);
        }
    }

    // initialize clusters to index of cavities
    std::vector<int> cluster_numbers;
    for (int i = 0; i < records.size(); i++) cluster_numbers.emplace_back(i);

    // buildClusters()

    for (int i=0; i<records.size(); i++) 
    for (int j=0; j<pairs.size(); j++) 
        if (pairs[j].first == i) cluster_numbers[pairs[j].second] = cluster_numbers[i];
/* was
    for (int i=0; i<number_of_cavities; i++) 
    for (int j=0; j<number_of_pairs; j++) 
        if (cavityA[j] == i) cluster_number[cavityB[j]] = cluster_number[i];
*/


    // build cluster objects 
    // cluster_numbers are the one-to-one mapping of cavity number to cluster number
    
    for (int i = 0; i < cluster_numbers.size(); i++) // for each cluster number
    {
        CavityCluster cluster; // create an empty cluster

        for (int j = 0; j < cluster_numbers.size(); j++) // find all cavities mapped to that cluster number
            if (cluster_numbers[j] == i) cluster.pushBack(records[i]);
        
        // size zero clusters get skipped and discarded
        if (cluster.getSize() > 0) clusters.push_back(cluster);
    }


    // deleteEmptyClusters()

/* was    
    int number_of_clusters=number_of_cavities;

    // get rid of empty clusters...
    // first loop over cluster numbers...
    for (int i=0; i<number_of_clusters;)
    {
        int cavs_in_cluster=0;
        for (int j=0; j<number_of_cavities; j++) if (cluster_number[j] == i) cavs_in_cluster++;

        if (cavs_in_cluster == 0)
        {
            // empty cluster, shift all higher cluster numbers down...
            for (j=0; j<number_of_cavities; j++) if (cluster_number[j] > i) cluster_number[j]--;
            number_of_clusters--;
        }
        else i++;
    }



    // sortClusters()

    double temp_x, temp_y, temp_z, temp_d, temp_cluster_no;
  
    for (int i=0; i<number_of_cavities-1; i++)
    for (int j=i; j<number_of_cavities; j++)
    {
        if (cluster_number[j] < cluster_number[i])
        {
            // swop the two cavities.

            vacuumms_float temp_x = x[i];
            vacuumms_float temp_y = y[i];
            vacuumms_float temp_z = z[i];
            vacuumms_float temp_d = d[i];
            vacuumms_float temp_cluster_no = cluster_number[i];

            x[i] = x[j];
            y[i] = y[j];
            z[i] = z[j];
            d[i] = d[j];
            cluster_number[i] = cluster_number[j];

            x[j] = temp_x;
            y[j] = temp_y;
            z[j] = temp_z;
            d[j] = temp_d;
            cluster_number[j] = temp_cluster_no;
        }
    }
*/

/* printCluster()
  int i, j;
  for (i=0; i<number_of_cavities; i++)
    printf("%d\t%05d\t%lf\t%lf\t%lf\t%lf\n", cluster_number[i], i, x[i], y[i], z[i], d[i]);
*/

    return clusters;

} // end generateClusters()


int CavityCluster::getSize()
{
    return configuration.getSize();
}


