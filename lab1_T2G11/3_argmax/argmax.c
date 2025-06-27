#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.14159265358979323846
#include <math.h>
#include <omp.h>

void initialize(double *v, int N) {
    for (int i = 0; i < N; i++) {
        v[i] = (1 - pow(0.5 - (double)i/(double)N, 2)) * cos(2*M_PI*100* (i - 0.5)/N);
    }
}

// computes the argmax sequentially with a for loop
void argmax_seq(double *v, int N, double *m, int *idx_m) {
    *m = v[0];
    *idx_m = 0;
    for (int i = 1; i < N; i++) {
        if (v[i] > *m) {
            *m = v[i];
            *idx_m = i;
        }
    }
}

// computes the argmax in parallel with a for loop
void argmax_par(double *v, int N, double *m, int *idx_m) {
    *m = v[0];
    *idx_m = 0;

    #pragma omp parallel
    {
        double thread_m = v[0];
        int thread_idx = 0;
        
        #pragma omp for
        for (int i = 1; i < N; i++) {
            if (v[i] > thread_m) {
                thread_m = v[i];
                thread_idx = i;
            }
        }

        #pragma omp critical
        {
            if (thread_m > *m) {
                *m = thread_m;
                *idx_m = thread_idx;
            }
        }
    }
}

// computes the argmax recursively and sequentially
void argmax_recursive(double *v, int N, double *m, int *idx_m, int K) {
    if(N<K){
        argmax_seq(v,N,m,idx_m);
        return;
    }
    else{
        double max1;
        double max2;
        int idx_1;
        int idx_2;
        int mid = N/2;

        argmax_recursive(v, mid, &max1, &idx_1, K);
        argmax_recursive(v + mid, N - mid, &max2, &idx_2, K);

        if(max1>max2){
            *m = max1;
            *idx_m = idx_1;
        }else{
            *m = max2;
            *idx_m = idx_2 + mid;
        }
    }

}

// computes the argmax recursively and in parallel using tasks
void argmax_recursive_tasks(double *v, int N, double *m, int *idx_m, int K) {
    if(N<=K){
        argmax_seq(v,N,m,idx_m);
        return;
    }
    else{
        double max1;
        double max2;
        int idx_1;
        int idx_2;
        int mid = N/2;
        
        #pragma omp task shared(max1,idx_1)
        argmax_recursive(v,mid,&max1,&idx_1,K);
        #pragma omp task shared(max2,idx_2)
        argmax_recursive(v+mid,N-mid,&max2,&idx_2,K);
        
        #pragma omp taskwait 
        if(max1>max2){
            *m = max1;
            *idx_m = idx_1;
        }else{
            *m = max2;
            *idx_m = idx_2 + mid;
        }
    }
}

int main(int argc, char* argv[]){
    int K;
    int nthreads;

    if(argc == 3) {
        K = atoi(argv[2]);
        nthreads = atoi(argv[1]);
    } else {
        return -1;
    }

    int N = pow(4096,2);
    
    double time;

    omp_set_num_threads(nthreads);

    double *v = (double *)malloc(N*sizeof(double));
    initialize(v,N);

    #pragma omp parallel
        #pragma omp single
        printf("Running argmax with K = %d using %d threads\n",K,omp_get_num_threads());


    ////////////////////////////////////////////////////////////////
    // Sequential
    ////////////////////////////////////////////////////////////////
    printf("Sequential for argmax: ");
    double seq_m;
    int seq_idx_m;
    
    time = omp_get_wtime();
    argmax_seq(v,N,&seq_m,&seq_idx_m);
    time = omp_get_wtime() - time;
    printf("m = %lf, idx_m = %d, time = %lf\n", seq_m, seq_idx_m,time);

    ////////////////////////////////////////////////////////////////
    // Parallel
    ////////////////////////////////////////////////////////////////
    printf("Parallel for argmax: ");
    double par_m;
    int par_idx_m;
    
    time = omp_get_wtime();
    argmax_par(v,N,&par_m,&par_idx_m);
    time = omp_get_wtime() - time;
    printf("m = %lf, idx_m = %d, time = %lf\n", par_m, par_idx_m,time);

    ////////////////////////////////////////////////////////////////
    // Sequential Recursive
    ////////////////////////////////////////////////////////////////
    printf("Sequential recursive argmax: ");
    double rec_m;
    int rec_idx_m;
    
    time = omp_get_wtime();
    argmax_recursive(v,N,&rec_m,&rec_idx_m,K);
    time = omp_get_wtime() - time;
    printf("m = %lf, idx_m = %d, time = %lf\n", rec_m, rec_idx_m,time);

    ////////////////////////////////////////////////////////////////
    // Parallel Recursive
    ////////////////////////////////////////////////////////////////
    printf("Parallel recursive argmax: ");
    double task_m;
    int task_idx_m;
    
    time = omp_get_wtime();
    #pragma omp parallel
        #pragma omp single
        argmax_recursive_tasks(v,N,&task_m,&task_idx_m,K);
    time = omp_get_wtime() - time;
    printf("m = %lf, idx_m = %d, time = %lf\n", task_m, task_idx_m,time);



    free(v);

    return 0;
}