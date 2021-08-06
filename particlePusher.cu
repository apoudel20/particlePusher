#include <stdio.h>
#include <iostream>

#include <cstdlib> // for atoi

#include <fstream>

#define GRID_CENTER_OFFSET .25
#define GRID_SIZE 100

__device__ float vectorFieldX(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_x = real_x/10;//pow(real_x,2) * real_y;
    return new_vel_x;
}
__device__ float vectorFieldY(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_y = 1/(exp(pow(0.1*real_x,4)))*1/(exp(pow(0.1*real_z,4)));//1/sqrt(1+pow((real_x),2) * pow((real_z),2));
    return new_vel_y;
}
__device__ float vectorFieldZ(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_z = real_z/10;//real_x * real_y * real_z;
    return new_vel_z;
}

__global__ 
void rayTrace(int N, float *real_x, float *real_y, float *real_z,float *vel_x, float *vel_y, float *vel_z, int timesteps){//, curandState *state){
    int rayIndex = blockIdx.x*blockDim.x + threadIdx.x;

    if(rayIndex < N){
        real_x[rayIndex] = ((float)rayIndex)/N*.5*GRID_CENTER_OFFSET*GRID_SIZE * cospi(60*((float)rayIndex)/N);
        real_y[rayIndex] = 0;
        real_z[rayIndex] = ((float)rayIndex)/N*.5*GRID_CENTER_OFFSET*GRID_SIZE * sinpi(60*((float)rayIndex)/N);
    }
    float d_vel_x, d_vel_y, d_vel_z;

    d_vel_x = vectorFieldX(real_x[rayIndex], real_y[rayIndex], real_z[rayIndex], vel_x[rayIndex],vel_y[rayIndex],vel_z[rayIndex],timesteps);
    d_vel_y = vectorFieldY(real_x[rayIndex], real_y[rayIndex], real_z[rayIndex], vel_x[rayIndex],vel_y[rayIndex],vel_z[rayIndex],timesteps);
    d_vel_z = vectorFieldZ(real_x[rayIndex], real_y[rayIndex], real_z[rayIndex], vel_x[rayIndex],vel_y[rayIndex],vel_z[rayIndex],timesteps);    

    float inter_vel_x, inter_vel_y, inter_vel_z; // intermediate values before calculating boundary conditions for reflection off walls

    inter_vel_x = vel_x[rayIndex] + d_vel_x + real_x[rayIndex];
    inter_vel_y = vel_y[rayIndex] + d_vel_y + real_y[rayIndex];
    inter_vel_z = vel_z[rayIndex] + d_vel_z + real_z[rayIndex];

    // boundary reflection conditions
    vel_x[rayIndex] += d_vel_x;
    if(inter_vel_x > GRID_SIZE/8){
        real_x[rayIndex] = 2 * GRID_SIZE/8 - inter_vel_x;
        vel_x[rayIndex] = -vel_x[rayIndex];
    }else if(inter_vel_x < -GRID_SIZE/8){
        real_x[rayIndex] = -GRID_SIZE/8 - inter_vel_x;
        vel_x[rayIndex] = -vel_x[rayIndex];
    }else{
        real_x[rayIndex] += vel_x[rayIndex];
    }

    vel_y[rayIndex] += d_vel_y;
    if(inter_vel_y > GRID_SIZE/10){
        real_y[rayIndex] = 2 * GRID_SIZE/10 - inter_vel_y;
        vel_y[rayIndex] = -vel_y[rayIndex];
    }else if(inter_vel_y < 0){
        real_y[rayIndex] = -inter_vel_y;
        vel_y[rayIndex] = -vel_y[rayIndex];
    }else{
        real_y[rayIndex] += vel_y[rayIndex];
    }

    vel_z[rayIndex] += d_vel_z;
    if(inter_vel_z > GRID_SIZE/8){
        real_z[rayIndex] = 2 * GRID_SIZE/8 - inter_vel_z;
        vel_z[rayIndex] = -vel_z[rayIndex];
    }else if(inter_vel_z < -GRID_SIZE/8){
        real_z[rayIndex] = -GRID_SIZE/8 - inter_vel_z;
        vel_z[rayIndex] = -vel_z[rayIndex];
    }else{
        real_z[rayIndex] += vel_z[rayIndex];
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

    
    float *real_x, *real_y, *real_z,*dev_real_x, *dev_real_y, *dev_real_z; // R3 coordinates to be mapped onto the gridspace 
    float *vel_x, *vel_y, *vel_z, *dev_vel_x, *dev_vel_y, *dev_vel_z; // current velocities for the particles

    // allocate space for position arrays    
    real_x = (float*)malloc(N * sizeof(float));
    real_y = (float*)malloc(N * sizeof(float));
    real_z = (float*)malloc(N * sizeof(float));
    
    // allocate space for velocity arrays
    vel_x = (float*)malloc(N * sizeof(float));
    vel_y = (float*)malloc(N * sizeof(float));
    vel_z = (float*)malloc(N * sizeof(float));

    // allocate device space for position arrays
    cudaMalloc(&dev_real_x, N * sizeof(float));
    cudaMalloc(&dev_real_y, N * sizeof(float));
    cudaMalloc(&dev_real_z, N * sizeof(float));

    // allocate device space for velocity arrays
    cudaMalloc(&dev_vel_x, N * sizeof(float));
    cudaMalloc(&dev_vel_y, N * sizeof(float));
    cudaMalloc(&dev_vel_z, N * sizeof(float));
    
    // initialize arrays to 0
    cudaMemset(dev_real_x, 0, N * sizeof(float));
    cudaMemset(dev_real_y, 0, N * sizeof(float));
    cudaMemset(dev_real_z, 0, N * sizeof(float));
    cudaMemset(dev_vel_x, 0, N * sizeof(float));
    cudaMemset(dev_vel_y, 0, N * sizeof(float));
    cudaMemset(dev_vel_z, 0, N * sizeof(float));

    for(int t = 0; t < timesteps; t++){
        rayTrace<<<(N+192)/256,256>>>(N, dev_real_x, dev_real_y, dev_real_z,dev_vel_x, dev_vel_y, dev_vel_z, t);//, devStates);
        cudaDeviceSynchronize();
    }

    // return results to host
    cudaMemcpy(real_x, dev_real_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_y, dev_real_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_z, dev_real_z, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_x, dev_vel_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_y, dev_vel_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_z, dev_vel_z, N * sizeof(float), cudaMemcpyDeviceToHost);


    std::ofstream outputWrite;
    outputWrite.open("frames/output_"+std::to_string(timesteps)+".csv");

    std::cout<<"Generating until timestep: "<<timesteps<<std::endl;
    for (int i = 0; i < N; i++) {
        outputWrite << real_x[i] << "," << real_y[i] << "," << real_z[i]<<std::endl;
    }

    outputWrite.close();

    cudaFree(dev_vel_z);
    cudaFree(dev_vel_y);
    cudaFree(dev_vel_x);    
    cudaFree(dev_real_z);
    cudaFree(dev_real_y);
    cudaFree(dev_real_x);

    free(vel_z);
    free(vel_y);
    free(vel_x);
    free(real_z);
    free(real_y);
    free(real_x);
    


}