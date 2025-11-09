#include "mpi.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct{ 
    uint64_t state;
    uint64_t inc;
} pcg32_random_t;

double pcg32_random(pcg32_random_t* rng){
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    uint32_t ran_int = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return (double)ran_int / (double)UINT32_MAX;
}

int main(int argc, char **argv){
    int rank,size;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N;
    long num_samples,SEED;
    
    if(argc==4){
        N = atoi(argv[1]);
        num_samples = atol(argv[2]);
        SEED = atol(argv[3]);
    }else{
        N = 3;
        num_samples = 1000000;
        SEED = time(NULL);
    }

    double start, end;
    double ratio;
    double x;
    start = MPI_Wtime();
       

    pcg32_random_t rng;
    rng.state = SEED + rank;
    rng.inc = (rank << 16) | 0x3039;
    ratio = pow(M_PI,N/2)*1/(tgamma((N/2)+1)*pow(2,N));

    long count = 0;
    long samples = num_samples/size;
    long residu = num_samples%size;

    if(rank<residu){
        samples+=1;
    }

    for(long i=0;i<samples;i++){
        double x2 = 0.0;
        for(int n=0; n<N ; n++){
            x = pcg32_random(&rng);
            double aux = 2*x - 1;
            x2 += aux*aux;
        }
        double norm = sqrt(x2);
        if(norm<=1.0){
            count++;
        } 
    }

    long global_count = 0;
    MPI_Reduce(&count,&global_count,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    end = MPI_Wtime()-start;

    double slowest = 0.0;
    MPI_Reduce(&end,&slowest,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(rank==0){
        double global_ratio = (double)global_count/num_samples;
        double error = fabs(global_ratio-ratio);
        printf("Monte Carlo sphere/cube ratio estimation\n");
        printf("N: %ld samples, d: %d, seed %d, size: %d\n",num_samples,N,SEED,size);
        printf("Ratio = %lf Err: %lf\n",global_ratio,error);
        printf("Elapsed time: %lf\n",slowest);
    }

    MPI_Finalize();
    
    return 0;

}

