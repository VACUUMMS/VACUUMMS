/* fvi.cu */

#include <cuda_runtime.h>
#include <vector>

#include <vacuumms/fvi.hh>

// This kernel will calculate attraction, repulsion, energy, and FVI,
// for any metric where the device pointer passed is not NULL.
__global__ void Kernel(
    ConfigurationRecord* d_configuration, 
    int     n_records, 
    vacuumms_float   box_x, 
    vacuumms_float   box_y, 
    vacuumms_float   box_z,
    size_t           dim_x,
    size_t           dim_y,
    size_t           dim_z,
    vacuumms_float *d_attraction,
    vacuumms_float *d_repulsion,
    vacuumms_float *d_energy,
    vacuumms_float *d_FVI)
{
    // blockIdx and blockDim values are provided by CUDA
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idz = blockIdx.z * blockDim.z + threadIdx.z;

    vacuumms_float repulsion=0;
    vacuumms_float attraction=0;
    vacuumms_float sigma_over_r_sq;
    vacuumms_float dx, dy, dz, dd;
    vacuumms_float f_resolution_x = box_x / dim_x;
    vacuumms_float f_resolution_y = box_y / dim_y;
    vacuumms_float f_resolution_z = box_z / dim_z;

    vacuumms_float cuda_x = idx * f_resolution_x;
    vacuumms_float cuda_y = idy * f_resolution_y;
    vacuumms_float cuda_z = idz * f_resolution_z;

    vacuumms_float sigma_probe = 0.0f;
    vacuumms_float epsilon_probe = 1.0f;

    // evaluate energy at (cuda_x, cuda_y, cuda_z);
    for (int i=0; i< n_records; i++) 
    {
        // Lorentz-Berthelot combining rules
        vacuumms_float sigma_ij = 0.5 * (d_configuration[i].sigma + sigma_probe);
        vacuumms_float sigma_ij_sq = sigma_ij * sigma_ij;
        vacuumms_float epsilon_ij = sqrt(d_configuration[i].sigma * epsilon_probe);

        // loop over mirror boxes
        for (int l=-1; l<=1; l++) 
        for (int m=-1; m<=1; m++) 
        for (int n=-1; n<=1; n++) 
        {
            // central atom
            dx = l * box_x + d_configuration[i].x - cuda_x;
            dy = m * box_y + d_configuration[i].y - cuda_y;
            dz = n * box_z + d_configuration[i].z - cuda_z;
            dd = dx*dx + dy*dy + dz*dz; 
   
            sigma_over_r_sq = sigma_ij_sq / dd; 
            vacuumms_float sigma_over_r_6 = sigma_over_r_sq * sigma_over_r_sq * sigma_over_r_sq;
            vacuumms_float sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
            repulsion += d_configuration[i].epsilon * sigma_over_r_12;
            attraction += d_configuration[i].epsilon * sigma_over_r_6;
        }
    } 

    size_t which = idx * dim_x * dim_y + idy * dim_y + idz;
    if (d_attraction != nullptr) d_attraction[which] = 4 * attraction;
    if (d_repulsion != nullptr) d_repulsion[which] = 4 * repulsion;
    if (d_energy != nullptr) d_energy[which] = 4 * repulsion - 4 * attraction;
    if (d_FVI != nullptr) d_FVI[which] = exp(-4 * repulsion);

} // end of Kernel


// This is the routine which is exposed in the API
void FVIX::executeMask(int mask)
{
    size_t n_records = c.getSize(); 
    size_t array_size = dimensions[0] * dimensions[1] * dimensions[2];

    std::vector<ConfigurationRecord> h_records;
    for (int i=0; i<n_records; i++) h_records.push_back(ConfigurationRecord(c.recordAt(i)));

    /* allocate for return values on device */

    vacuumms_float* d_repulsion = NULL;
    if (mask & FVIX_REPULSION)
        for(cudaError_t err = cudaErrorUnknown; 
            err != cudaSuccess; 
            err = cudaMalloc( &d_repulsion, array_size * sizeof(vacuumms_float)));

    vacuumms_float* d_attraction = NULL;
    if (mask & FVIX_ATTRACTION)
        for(cudaError_t err = cudaErrorUnknown; 
            err != cudaSuccess; 
            err = cudaMalloc( &d_attraction, array_size * sizeof(vacuumms_float)));

    vacuumms_float* d_energy = NULL;
    if (mask & FVIX_ENERGY)
        for(cudaError_t err = cudaErrorUnknown; 
            err != cudaSuccess; 
            err = cudaMalloc( &d_energy, array_size * sizeof(vacuumms_float)));

    vacuumms_float* d_FVI = NULL;
    if (mask & FVIX_FVI)
        for(cudaError_t err = cudaErrorUnknown; 
            err != cudaSuccess; 
            err = cudaMalloc( &d_FVI, array_size * sizeof(vacuumms_float)));

    /* malloc, copy, and sync config records */

    ConfigurationRecord *d_records;

    for(cudaError_t err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMalloc( &d_records, sizeof(ConfigurationRecord) * n_records));

    for(cudaError_t err = cudaErrorUnknown; 
        err != cudaSuccess; 
        err = cudaMemcpy( d_records, h_records.data(), h_records.size() * sizeof(ConfigurationRecord), cudaMemcpyHostToDevice ));

    cudaDeviceSynchronize(); // block until the device has completed
    cudaError_t last = cudaGetLastError();
    if (last != cudaSuccess) printf("%s\n", cudaGetErrorString(last)); 


    dim3 dimBlock(8, 8, 8);
    dim3 dimGrid(
        (dimensions[0] + dimBlock.x -1) / dimBlock.x,
        (dimensions[1] + dimBlock.y -1) / dimBlock.y,
        (dimensions[2] + dimBlock.z -1) / dimBlock.z
    );

    if ( (dimensions[0] < 8) || (dimensions[1] < 8) || (dimensions[2] < 8) )
    {
        printf("FVIX::Kernel: need all dimensions to be >= 8\n");
        return;
    }

    Kernel<<< dimGrid, dimBlock >>>(
        d_records, 
        n_records, 
        c.box_x, 
        c.box_y, 
        c.box_z, 
        dimensions[0], 
        dimensions[1], 
        dimensions[2], 
        d_attraction, 
        d_repulsion, 
        d_energy, 
        d_FVI);

    cudaDeviceSynchronize(); // block until the device has completed
    last = cudaGetLastError();
    if (last != cudaSuccess) printf("%s\n", cudaGetErrorString(last)); 

/*
    // retrieve result
    attraction.resize(array_size);
    repulsion.resize(array_size);
    energy.resize(array_size);
    FVI.resize(array_size);
*/

    cudaError_t err;

    if (mask & FVIX_ATTRACTION)
    {
        attraction.resize(array_size);
        err = cudaMemcpy(attraction.data(), d_attraction, array_size * sizeof(vacuumms_float), cudaMemcpyDeviceToHost );
        if (err != cudaSuccess) std::cerr <<  "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
        cudaDeviceSynchronize(); // block until the device has completed
        err = cudaGetLastError();
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 
    }
      
    if (mask & FVIX_REPULSION)
    {
        repulsion.resize(array_size);
        err = cudaMemcpy(repulsion.data(), d_repulsion, array_size * sizeof(vacuumms_float), cudaMemcpyDeviceToHost );
        if (err != cudaSuccess) std::cerr <<  "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
        cudaDeviceSynchronize(); // block until the device has completed
        err = cudaGetLastError();
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 
    }

      
//if (d_energy != NULL)
    if (mask & FVIX_ENERGY)
    {
        energy.resize(array_size);
        err = cudaMemcpy(energy.data(), d_energy, array_size * sizeof(vacuumms_float), cudaMemcpyDeviceToHost );
        if (err != cudaSuccess) std::cerr <<  "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
        cudaDeviceSynchronize(); // block until the device has completed
        err = cudaGetLastError();
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 
    }

      
    if (mask & FVIX_FVI)
    {
        FVI.resize(array_size);
        err = cudaMemcpy(FVI.data(), d_FVI, array_size * sizeof(vacuumms_float), cudaMemcpyDeviceToHost );
        if (err != cudaSuccess) std::cerr <<  "cudaMemcpy failed: " << cudaGetErrorString(err) << std::endl;
        cudaDeviceSynchronize(); // block until the device has completed
        err = cudaGetLastError();
        if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 
    }


    cudaFree(d_records);
    cudaFree(d_attraction);
    cudaFree(d_repulsion);
    cudaFree(d_energy);
    cudaFree(d_FVI);

    cudaDeviceSynchronize(); // block until the device has completed
    err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err)); 

}
