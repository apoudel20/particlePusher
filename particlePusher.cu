#include <stdio.h>
#include <iostream>

#include <cstdlib> // for atoi

#include <fstream>

#define GRID_CENTER_OFFSET .25
#define GRID_SIZE 100



__global__ 
void rayTrace(int N, float *real_x, float *real_y, float *real_z, float timesteps){//, curandState *state){
    int rayIndex = blockIdx.x*blockDim.x + threadIdx.x;

    if(rayIndex < N){
        real_x[rayIndex] = ((float)rayIndex)/N*.5*GRID_CENTER_OFFSET*GRID_SIZE * cospi(60*((float)rayIndex)/N);
        real_y[rayIndex] = 0;
        real_z[rayIndex] = ((float)rayIndex)/N*.5*GRID_CENTER_OFFSET*GRID_SIZE * sinpi(60*((float)rayIndex)/N);
    }

    real_y[rayIndex] += timesteps;
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

    rayTrace<<<(N+192)/256,256>>>(N, dev_real_x, dev_real_y, dev_real_z, timesteps);//, devStates);

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

    cudaFree(dev_real_z);
    cudaFree(dev_real_y);
    cudaFree(dev_real_x);


}