#include <stdio.h>
#include <iostream>

#include <cstdlib> // for atoi

#include <fstream>

#include <curand_kernel.h>

#define GRID_CENTER_OFFSET .25
#define GRID_SIZE 100

// __global__
// void saxpy(int n, float a, float *x, float *y)
// {
//     int i = blockIdx.x*blockDim.x + threadIdx.x;
//     if (i < n) y[i] = a*x[i] + y[i];
// }

// __device__ 
// void initializeParticlePosition(int n, float *real_x, float *real_y, float *real_z){


// }

__global__ void setup_kernel(curandState *state)
{
    int id =  blockIdx.x*blockDim.x + threadIdx.x;
    /* Each thread gets same seed , a different sequence
    number , no offset */
    curand_init (1234 , id, 0, &state[id]);
}

__global__ void generate_kernel(curandState *state ,int *result){
    int id = threadIdx.x + blockIdx.x * 64;
    int count = 0;
    unsigned int x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[id];
    /* Generate pseudo -random unsigned ints */
    for(int n = 0; n < 100000; n++) {
        x = curand (& localState);
        /* Check if low bit set */
        if(x & 1) {
            count ++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ 
void rayTrace(int N, float *real_x, float *real_y, float *real_z, float timesteps, curandState *state){
    int rayIndex = blockIdx.x*blockDim.x + threadIdx.x;
    curandState localState = state[rayIndex];
    float random = curand_uniform(&localState);
    printf("%f \n--\n",random);
    if(rayIndex < N){
        real_x[rayIndex] = (GRID_SIZE/2 + (GRID_CENTER_OFFSET * GRID_SIZE * cospi(20 * ((float)rayIndex)/N)))* random*10000;
        real_y[rayIndex] = 0;
        real_z[rayIndex] = (GRID_SIZE/2 + (GRID_CENTER_OFFSET * GRID_SIZE * sinpi(20 * ((float)rayIndex)/N)))* random*10000;
    }
    // state[rayIndex] = localState;

    for(int t = 0; t < timesteps; t++){
        real_y[rayIndex] += 1;
    }
}

int main(int argc, char *argv[]) 
{
    const int N = pow(GRID_SIZE, 3);

    int timesteps = 0;

    if(argc == 2){
        timesteps = atoi(argv[1]);
    }else{
        timesteps = 0;
    }

    
    float *real_x, *real_y, *real_z,*dev_real_x, *dev_real_y, *dev_real_z;
    
    real_x = (float*)malloc(N * sizeof(float));
    real_y = (float*)malloc(N * sizeof(float));
    real_z = (float*)malloc(N * sizeof(float));
    
    cudaMalloc(&dev_real_x, N * sizeof(float));
    cudaMalloc(&dev_real_y, N * sizeof(float));
    cudaMalloc(&dev_real_z, N * sizeof(float));
    
    cudaMemset(dev_real_x, 0, N * sizeof(float));
    cudaMemset(dev_real_y, 0, N * sizeof(float));
    cudaMemset(dev_real_z, 0, N * sizeof(float));
    // cuRAND setup
    curandState *devStates;
    cudaMalloc(&devStates, 2 * N * sizeof(float));

    setup_kernel<<<(2*N+128)/256,256>>>(devStates);
    cudaDeviceSynchronize();

    rayTrace<<<(N+192)/256,256>>>(N, dev_real_x, dev_real_y, dev_real_z, timesteps, devStates);

    cudaMemcpy(real_x, dev_real_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_y, dev_real_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_z, dev_real_z, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    std::ofstream outputWrite;
    outputWrite.open("frames/output_"+std::to_string(timesteps)+".csv");

    std::cout<<"Generating until timestep: "<<timesteps<<std::endl;
    for (int i = 0; i < N; i++) {
        outputWrite << real_x[i] << "," << real_y[i] << "," << real_z[i]<<std::endl;
    }

    outputWrite.close();

    cudaFree(devStates);
    cudaFree(dev_real_z);
    cudaFree(dev_real_y);
    cudaFree(dev_real_x);


}