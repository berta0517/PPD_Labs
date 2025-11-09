#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#define BLOCK_SIZE 128
#define CUDA_CHECK(call)                                                \
    do {                                                                    \
        cudaError_t err = (call);                                           \
        if (err != cudaSuccess) {                                           \
            fprintf(stderr, "CUDA error at %s:%d: %s\n",                    \
                    __FILE__, __LINE__, cudaGetErrorString(err));          \
            return 1;                                                       \
        }                                                                   \
    } while (0)

// Kernel de suma de vectores
__global__ void vecadd_cuda(int *a, int *b, int *c, int n) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < n)
        c[index] = a[index] + b[index];
}



// Guarda la salida en un archivo (puedes personalizar esto)
void save_output(int *c, int n) {
    FILE *f = fopen("output.txt", "w");
    for (int i = 0; i < n; i++) {
        fprintf(f, "%d\n", c[i]);
    }
    fclose(f);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <vector_size>\n", argv[0]);
        return 1;
    }
    float copy_time, kernel_time, device_host_time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    int N = atoi(argv[1]);
    int size = N * sizeof(int);

    // Reserva de memoria en el host
    int *a = (int *)malloc(size);
    int *b = (int *)malloc(size);
    int *c = (int *)malloc(size);

    // Inicialización
    for (int i = 0; i < N; i++)
    {
        a[i] = i;
        b[i] = 2*(N-i);
    }
    // Reserva de memoria en el device
    int *d_a, *d_b, *d_c;
    CUDA_CHECK(cudaMalloc((void **)&d_a, size));
    CUDA_CHECK(cudaMalloc((void **)&d_b, size));
    CUDA_CHECK(cudaMalloc((void **)&d_c, size));

    // Copiar datos al dispositivo
    CUDA_CHECK(cudaEventRecord(start));
    CUDA_CHECK(cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&copy_time, start, stop));
    printf("Time to copy data to device: %f ms\n", copy_time);

    // Configurar y lanzar el kernel
    int number_of_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CUDA_CHECK(cudaEventRecord(start));
    vecadd_cuda<<<number_of_blocks, BLOCK_SIZE>>>(d_a, d_b, d_c, N);
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&kernel_time, start, stop));
    printf("Kernel execution time: %f ms\n", kernel_time);

    // Esperar a que termine y copiar el resultado al host
    CUDA_CHECK(cudaEventRecord(start));
    CUDA_CHECK(cudaMemcpy(c, d_c, size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&device_host_time, start, stop));
    printf("Time to copy data back to host: %f ms\n", device_host_time);

    int valid = 1;
    for (int i = 0; i < N; i++)
    {
        if (fabs(c[i] - (2*N-i)) > 10e-6){
            valid = 0;
            break;
        }
    }
    if(valid){
        printf("Validation successful\n");
    }
    else{
        printf("Validation failed\n");
    }

    // Liberar memoria
    cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
    free(a); free(b); free(c);

    return 0;
}
