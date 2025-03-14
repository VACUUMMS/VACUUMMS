/* fvi.cu */

#include <cuda_runtime.h>
#include <vector>
#include <vacuumms/cuda.h>
#include <vacuumms/cuda.hh>

// This is the kernel, called by the GFGToRepulsionX() functions, 
// which are, in turn, exposed as the API.
__global__ void EnergyKernel16_612(
    ConfigurationRecord* d_configuration, 
    int     n_records, 
    float   box_x, 
    float   box_y, 
    float   box_z,
    vacuumms_EnergyArray16 *d_repulsion)
{
    // blockIdx values are provided by CUDA
    unsigned int idx = blockIdx.x;
    unsigned int idy = blockIdx.y;
    unsigned int idz = threadIdx.x;

    float repulsion=0;
//    float attraction=0;
    float sigma_over_r_sq;
    float dx, dy, dz, dd;
    float f_resolution_x = box_x / 16.0f;
    float f_resolution_y = box_y / 16.0f;
    float f_resolution_z = box_z / 16.0f;

    float cuda_x = idx * f_resolution_x;
    float cuda_y = idy * f_resolution_y;
    float cuda_z = idz * f_resolution_z;

    // evaluate energy at (cuda_x, cuda_y, cuda_z);
    for (int i=0; i< n_records; i++) 
    {
        // central atom
        dx = d_configuration[i].x - cuda_x;
        dy = d_configuration[i].y - cuda_y;
        dz = d_configuration[i].z - cuda_z;
        dd = dx*dx + dy*dy + dz*dz; 
        sigma_over_r_sq = d_configuration[i].sigma 
                        * d_configuration[i].sigma 
                        / dd; 
        repulsion += d_configuration[i].epsilon 
                  * sigma_over_r_sq 
                  * sigma_over_r_sq 
                  * sigma_over_r_sq 
                  * sigma_over_r_sq 
                  * sigma_over_r_sq 
                  * sigma_over_r_sq;

 //       attraction += d_configuration->atom[i].epsilon 
 //                 * sigma_over_r_sq 
 //                 * sigma_over_r_sq 
 //                 * sigma_over_r_sq;
    } 

// If NULL pointers are passed for the attraction or repulsion, no values are returned.
//    if (d_attraction) d_attraction->energy[idx][idy][idz] = 4 * attraction;
//    if (d_repulsion) d_repulsion->energy[idx][idy][idz] = 4 * repulsion;
//    if (d_total) d_total->energy[idx][idy][idz] = 4 * repulsion - 4 * attraction;
    d_repulsion->energy[idx][idy][idz] = 4 * repulsion;
}



vacuumms_EnergyArray16* calculateRepulsions(Configuration gfg)
{
    vacuumms_EnergyArray16 	*d_repulsion;

    int     n_records=gfg.getSize(); 

    // replicate the gfg. FTW need to come back to this, will have overly empty edges without
    // vacuumms_GFG65536 *h_configuration = replicateGFG65536(gfg); 



/*
    // and cross-parameterize use 6-12 rule
    for (int n=0; n<gfg->n_atoms; n++)
    {
        h_configuration->atom[n].sigma = 0.5f * (sigma + h_configuration->atom[n].sigma);
        h_configuration->atom[n].epsilon = sqrt(epsilon * h_configuration->atom[n].epsilon);
    }
*/

    // Push atoms into the container to be passed. Can add replication here later.
    std::vector<ConfigurationRecord> h_records;
    for (int i=0; i<gfg.getSize(); i++)
    {
        h_records.push_back(ConfigurationRecord(gfg.recordAt(i)));
printf("h_records = %d\n", h_records.size());
    }

printf("box = %f, %f, %f\n", gfg.box_x, gfg.box_y, gfg.box_z);

    cudaError_t err;
    /* allocate for energy array and configuration on device */
    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
//        err = cudaMalloc( (void **) &d_repulsion, sizeof(vacuumms_EnergyArray16)));
        err = cudaMalloc( &d_repulsion, sizeof(vacuumms_EnergyArray16)));
printf("successfully allocated d_repulsion\n");

    ConfigurationRecord *d_records;

    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
//        err = cudaMalloc( (void **) &d_configuration, sizeof(ConfigurationRecord) * n_records));
        err = cudaMalloc( &d_records, sizeof(ConfigurationRecord) * n_records));
printf("successfully allocated d_records\n");

    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
//        err = cudaMemcpy( h_records.data(), d_records, sizeof(ConfigurationRecord) * n_records, cudaMemcpyHostToDevice ));
        err = cudaMemcpy( d_records, h_records.data(), h_records.size() * sizeof(ConfigurationRecord), cudaMemcpyHostToDevice ));
printf("successfully copied h_records/d_records\n");

printf("sync device\n");

    cudaDeviceSynchronize(); // block until the device has completed
    err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 

    dim3 dimGrid(16, 16);
    dim3 dimBlock(16, 1, 1);

printf("running kernel\n");

    //EnergyKernel256_612<<< dimGrid, dimBlock >>>(d_configuration, NULL, d_repulsion, NULL);
    EnergyKernel16_612<<< dimGrid, dimBlock >>>(d_records, n_records, gfg.box_x, gfg.box_y, gfg.box_z, d_repulsion);

printf("sync device\n");

    cudaDeviceSynchronize(); // block until the device has completed
    err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 

printf("retrieve results\n");

    // retrieve result
    vacuumms_EnergyArray16 *h_repulsion = (vacuumms_EnergyArray16 *)malloc(sizeof(vacuumms_EnergyArray16));
    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMemcpy(h_repulsion, d_repulsion, sizeof(vacuumms_EnergyArray16), cudaMemcpyDeviceToHost ));

printf("free device mem\n");

  // free device memory
    //cudaFree(d_configuration);
    cudaFree(d_records);
    cudaFree(d_repulsion);

// This is an object now and will go out of scope    
//    free(h_configuration); // free host memory for replicated configuration

    return h_repulsion;
}

