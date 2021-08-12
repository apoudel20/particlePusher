#include <stdio.h>
#include <iostream>

#include <cstdlib> // for atoi

#include <fstream>

#define GRID_CENTER_OFFSET .25
#define GRID_SIZE 100
typedef struct
{
    double* x;
    double* y;
    double* z;
}TimeLog;
TimeLog* new_cuda_TimeLog(int size)
{
    TimeLog* result = (TimeLog*)malloc(sizeof(TimeLog));
    printf("%p %p %p\n", &(result->x), &(result->y), &(result->z));
    cudaMalloc(&(result->x), sizeof(double)*size);
    cudaMalloc(&(result->y), sizeof(double)*size);
    cudaMalloc(&(result->z), sizeof(double)*size);
    return result;
}
void copy_from_device(TimeLog* target, double* host_x,
                                       double* host_y,
                                       double* host_z,
                                       int size)
{
    if(target == NULL)
    {
        return;
    }
    cudaMemcpy(host_x, target->x, sizeof(double)*size, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_y, target->y, sizeof(double)*size, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_z, target->z, sizeof(double)*size, cudaMemcpyDeviceToHost);
}
__device__ float vectorFieldX(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_x = real_x/10;//pow(real_x,2) * real_y;
    return new_vel_x;
}
__device__ float vectorFieldY(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_y = 1/(exp(pow(0.1*real_x,4)))*1/(exp(pow(0.1*real_z,3)));//1/sqrt(1+pow((real_x),2) * pow((real_z),2));
    return new_vel_y;
}
__device__ float vectorFieldZ(float real_x, float real_y, float real_z,float vel_x, float vel_y, float vel_z,int timesteps){
    float new_vel_z = real_z/10;//real_x * real_y * real_z;
    return new_vel_z;
}
#define LOG 1
__global__ 
void rayTrace(int N, float *real_x, float *real_y, float *real_z,float *vel_x, float *vel_y, float *vel_z, int step, int timesteps, double dt, TimeLog* logger){//, curandState *state){
    int rayIndex = blockIdx.x*blockDim.x + threadIdx.x;
    if(rayIndex >= N)
    {
        return;
    }
    double *logx, *logy, *logz;
    if(LOG)
    {
        int offset = timesteps*rayIndex + step;
        logx = logger->x + offset;
        logy = logger->y + offset;
        logz = logger->z + offset;
    }

    double this_x, this_y, this_z;
    this_x = real_x[rayIndex];
    this_y = real_y[rayIndex];
    this_z = real_z[rayIndex];

    float d_vel_x, d_vel_y, d_vel_z;
    double this_vx, this_vy, this_vz; // intermediate values before calculating boundary conditions for reflection off walls

    this_vx = vel_x[rayIndex];
    this_vy = vel_y[rayIndex];
    this_vz = vel_z[rayIndex];
    d_vel_x = vectorFieldX(this_x, this_y, this_z, this_vx, this_vy, this_vz, step);
    d_vel_y = vectorFieldY(this_x, this_y, this_z, this_vx, this_vy, this_vz, step);
    d_vel_z = vectorFieldZ(this_x, this_y, this_z, this_vx, this_vy, this_vz, step); 
    this_vx = this_vx + d_vel_x;// + real_x[rayIndex];
    this_vy = this_vy + d_vel_y;// + real_y[rayIndex];
    this_vz = this_vz + d_vel_z;// + real_z[rayIndex];
    double future_x, future_y, future_z;
    future_x = this_vx*dt + this_x;
    future_y = this_vy*dt + this_y;
    future_z = this_vz*dt + this_z;
    // boundary reflection conditions
    float boundx = GRID_SIZE*.125;
    float boundy = GRID_SIZE*.08;
    float boundz = GRID_SIZE*.125;
    int condx = (future_x > boundx || future_x < -1*boundx);
    int condy = (future_y > boundy || future_y < -1*boundy);
    int condz = (future_z > boundz || future_z < -1*boundz);
    this_vx = -1*this_vx*condx + this_vx*!condx;
    this_vy = -1*this_vy*condy + this_vy*!condy;
    this_vz = -1*this_vz*condz + this_vz*!condz;
    this_x += this_vx*dt;
    this_y += this_vy*dt;
    this_z += this_vz*dt;
    if(LOG)
    {
        *logx = this_x;
        *logy = this_y;
        *logz = this_z;
    }
    real_x[rayIndex] = this_x;
    real_y[rayIndex] = this_y;
    real_z[rayIndex] = this_z;
    vel_x[rayIndex] = this_vx;
    vel_y[rayIndex] = this_vy;
    vel_z[rayIndex] = this_vz;
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
    for(int i = 0; i < N; i++)
    {
        real_x[i] = ((float)i)/N*.35*GRID_CENTER_OFFSET*GRID_SIZE*cospi(.8*GRID_SIZE*((float)i)/N);
        real_y[i] = 1;
        real_z[i] = ((float)i)/N*.35*GRID_CENTER_OFFSET*GRID_SIZE*sinpi(.8*GRID_SIZE*((float)i)/N);

    }
    // allocate device space for position arrays
    cudaMalloc(&dev_real_x, N * sizeof(float));
    cudaMalloc(&dev_real_y, N * sizeof(float));
    cudaMalloc(&dev_real_z, N * sizeof(float));

    // allocate device space for velocity arrays
    cudaMalloc(&dev_vel_x, N * sizeof(float));
    cudaMalloc(&dev_vel_y, N * sizeof(float));
    cudaMalloc(&dev_vel_z, N * sizeof(float));
    
    // initialize arrays to 0
    cudaMemcpy(dev_real_x, real_x, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_real_y, real_y, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_real_z, real_z, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemset(dev_vel_x, 0, N * sizeof(float));
    cudaMemset(dev_vel_y, 0, N * sizeof(float));
    cudaMemset(dev_vel_z, 0, N * sizeof(float));
    TimeLog* logger = new_cuda_TimeLog(timesteps*N), *cu_logger;
    cudaMalloc(&cu_logger, sizeof(TimeLog));
    cudaMemcpy(cu_logger, logger, sizeof(TimeLog), cudaMemcpyHostToDevice);
    for(int i = 0; i < timesteps; i++)
    {
        rayTrace<<<((int)pow(GRID_SIZE, 3) + ((256 - (int)pow(GRID_SIZE, 3) % 256) % 256))/256,256>>>(N, dev_real_x, dev_real_y, dev_real_z,dev_vel_x, dev_vel_y, dev_vel_z, i, timesteps, 1e-1, cu_logger);//, devStates);
        cudaDeviceSynchronize();
    }
    double* rx = (double*)malloc(sizeof(double)*timesteps*N);
    double* ry = (double*)malloc(sizeof(double)*timesteps*N);
    double* rz = (double*)malloc(sizeof(double)*timesteps*N);
    copy_from_device(logger, rx, ry, rz, timesteps*N);

    // return results to host
    cudaMemcpy(real_x, dev_real_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_y, dev_real_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(real_z, dev_real_z, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_x, dev_vel_x, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_y, dev_vel_y, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_z, dev_vel_z, N * sizeof(float), cudaMemcpyDeviceToHost);

    std::ofstream outputWrite;
    std::cout<<"Generating until timestep: "<<timesteps<<std::endl;
    for(int i = 0; i < timesteps; i++)
    {
        outputWrite.open("frames/output_"+std::to_string(i)+".csv");

        for (int j = 0; j < N; j++) {
            outputWrite << rx[j*timesteps+i] << "," << ry[j*timesteps+i] << "," << rz[j*timesteps+i]<<std::endl;
        }
        outputWrite.close();
        printf("frames/output_%s.csv\n", std::to_string(i).c_str());
    }


    cudaFree(dev_vel_z);
    cudaFree(dev_vel_y);
    cudaFree(dev_vel_x);    
    cudaFree(dev_real_z);
    cudaFree(dev_real_y);
    cudaFree(dev_real_x);
    //cudaFree(cu_logger);
    //free(logger);
    free(rx);
    free(ry);
    free(rz);
    free(vel_x);
    free(vel_y);
    free(vel_z);
    free(real_z);
    free(real_y);
    free(real_x);
    


}
