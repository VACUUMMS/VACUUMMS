/* fvi.cu */

#include <cuda_runtime.h>
#include <vector>

#include <vacuumms/cuda.h>
#include <vacuumms/fvi.hh>

// This is the kernel, called by the GFGToRepulsionX() functions, 
// which are, in turn, exposed as the API.
__global__ void EnergyKernel16_612(
    ConfigurationRecord* d_configuration, 
    int     n_records, 
    vacuumms_float   box_x, 
    vacuumms_float   box_y, 
    vacuumms_float   box_z,
    vacuumms_EnergyArray16 *d_attraction,
    vacuumms_EnergyArray16 *d_repulsion,
    vacuumms_EnergyArray16 *d_total)
{
    // blockIdx values are provided by CUDA
    unsigned int idx = blockIdx.x;
    unsigned int idy = blockIdx.y;
    unsigned int idz = threadIdx.x;

    vacuumms_float repulsion=0;
    vacuumms_float attraction=0;
    vacuumms_float sigma_over_r_sq;
    vacuumms_float dx, dy, dz, dd;
    vacuumms_float f_resolution_x = box_x / 16.0f;
    vacuumms_float f_resolution_y = box_y / 16.0f;
    vacuumms_float f_resolution_z = box_z / 16.0f;

    vacuumms_float cuda_x = idx * f_resolution_x;
    vacuumms_float cuda_y = idy * f_resolution_y;
    vacuumms_float cuda_z = idz * f_resolution_z;

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
        vacuumms_float sigma_over_r_6 = sigma_over_r_sq * sigma_over_r_sq * sigma_over_r_sq;
        vacuumms_float sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
        repulsion += d_configuration[i].epsilon * sigma_over_r_12;
        attraction += d_configuration[i].epsilon * sigma_over_r_6;
    } 

    // If NULL pointers are passed for the attraction or repulsion, no values are returned.
    if (d_attraction) d_attraction->energy[idx][idy][idz] = 4 * attraction;
    if (d_repulsion) d_repulsion->energy[idx][idy][idz] = 4 * repulsion;
    if (d_total) d_total->energy[idx][idy][idz] = 4 * repulsion - 4 * attraction;
}


template <size_t resolution>
FVIArray<resolution>* calculateFVI(Configuration gfg)
{
    fprintf(stderr, "calculateFVI not implemented and resolution = %d.\n", resolution);
    return nullptr;
}


vacuumms_EnergyArray16* FVIX::calculateRepulsions(Configuration gfg)
{
    vacuumms_EnergyArray16 	*d_repulsion;

    int     n_records=gfg.getSize(); 

    // Push atoms into the container to be passed. Can add replication here later.
    std::vector<ConfigurationRecord> h_records;
    for (int i=0; i<gfg.getSize(); i++)
    {
        h_records.push_back(ConfigurationRecord(gfg.recordAt(i)));
    }

    cudaError_t err;
    /* allocate for energy array and configuration on device */
    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMalloc( &d_repulsion, sizeof(vacuumms_EnergyArray16)));

    ConfigurationRecord *d_records;

    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMalloc( &d_records, sizeof(ConfigurationRecord) * n_records));
printf("successfully allocated d_records\n");

    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMemcpy( d_records, h_records.data(), h_records.size() * sizeof(ConfigurationRecord), cudaMemcpyHostToDevice ));

    cudaDeviceSynchronize(); // block until the device has completed
    err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 

    dim3 dimGrid(16, 16);
    dim3 dimBlock(16, 1, 1);

    EnergyKernel16_612<<< dimGrid, dimBlock >>>(d_records, n_records, gfg.box_x, gfg.box_y, gfg.box_z, NULL, d_repulsion, NULL);


    cudaDeviceSynchronize(); // block until the device has completed
    err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 

    // retrieve result
    vacuumms_EnergyArray16 *h_repulsion = (vacuumms_EnergyArray16 *)malloc(sizeof(vacuumms_EnergyArray16));
    for(err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMemcpy(h_repulsion, d_repulsion, sizeof(vacuumms_EnergyArray16), cudaMemcpyDeviceToHost ));

    cudaFree(d_records);
    cudaFree(d_repulsion);

    return h_repulsion;
}

