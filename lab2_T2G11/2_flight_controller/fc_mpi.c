#include <stdio.h>

#include <assert.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "auxiliar.h"

/// TODO
/// Reading the planes from a file for MPI
void read_planes_mpi(const char* filename, PlaneList* planes, int* N, int* M, double* x_max, double* y_max, int rank, int size, int* tile_displacements)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    char line[MAX_LINE_LENGTH];
    int num_planes = 0;

    // Reading header
    fgets(line, sizeof(line), file);
    fgets(line, sizeof(line), file);
    sscanf(line, "# Map: %lf, %lf : %d %d", x_max, y_max, N, M);
    fgets(line, sizeof(line), file);
    sscanf(line, "# Number of Planes: %d", &num_planes);
    fgets(line, sizeof(line), file);

    for(int i = 0;i<=size;i++){
        tile_displacements[i] = i*((*M)*(*N))/size;
    }

    // Reading plane data
    int planes_read = 0;
    while (fgets(line, sizeof(line), file)) {
        int idx;
        double x, y, vx, vy;
        if (sscanf(line, "%d %lf %lf %lf %lf", &idx, &x, &y, &vx, &vy) == 5) {
            int index_i = get_index_i(x, *x_max, *N);
            int index_j = get_index_j(y, *y_max, *M);
            int map_index = get_index(index_i, index_j, *N, *M);

            int plane_rank = get_rank_from_index(map_index,tile_displacements,size);
            if(plane_rank == rank){
                insert_plane(planes, idx, map_index, rank, x, y, vx, vy);
                planes_read++;
            }
        }
    }
    fclose(file);
    
    int total_read = 0;
    MPI_Reduce(&planes_read,&total_read,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0){
        printf("Total planes read: %d\n",total_read);
    }
}

/// TODO
/// Communicate planes using mainly Send/Recv calls with default data types
void communicate_planes_send(PlaneList* list, int N, int M, double x_max, double y_max, int rank, int size, int* tile_displacements)
{

    int *send_count = calloc(size,sizeof(int));

    PlaneNode* current = list->head;
    while (current != NULL) {
        int index_i = get_index_i(current->x, x_max, N);
        int index_j = get_index_j(current->y, y_max, M);
        int plane_rank = get_rank_from_indices(index_i, index_j,  N,  M, tile_displacements, size);
        if(plane_rank != rank){
            send_count[plane_rank]++;
        }
        current = current->next;

    }

    double** send_buffers = malloc(size*sizeof(double*));
    for (int i=0;i<size;i++) {
        if (send_count[i] > 0) {
            send_buffers[i] = malloc(5* send_count[i] * sizeof(double));
        }else{
            send_buffers[i] = NULL;
        }
    }

    current = list->head;
    int* counters = calloc(size,sizeof(int));
    while (current != NULL) {
        int i = get_index_i(current->x, x_max, N);
        int j = get_index_j(current->y, y_max, M);
        int map_rank = get_rank_from_indices(i, j, N, M, tile_displacements, size);

        if (map_rank != rank) {
            int index = counters[map_rank]*5;
            send_buffers[map_rank][index] = (double)current->index_plane;
            send_buffers[map_rank][index+1] = current->x;
            send_buffers[map_rank][index+2] = current->y;
            send_buffers[map_rank][index+3] = current->vx;
            send_buffers[map_rank][index+4] = current->vy;

            counters[map_rank]++;
        }

        current = current->next;
    }

    int* recv_counts = malloc(size*sizeof(int));
    MPI_Alltoall(send_count,1,MPI_INT,recv_counts,1,MPI_INT, MPI_COMM_WORLD);
    
    MPI_Request* send_req = malloc(size * sizeof(MPI_Request));
    for (int i=0;i<size;i++) {
        if(send_count[i] > 0){
            MPI_Isend(send_buffers[i], 5*send_count[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req[i]);
        }else{
            send_req[i] = MPI_REQUEST_NULL;
        }
    }


    double** recv_buffers = malloc(size*sizeof(double*));
    for (int i=0;i<size;i++) {
        if (recv_counts[i] > 0) {
            recv_buffers[i] = malloc(5*recv_counts[i]*sizeof(double));
            MPI_Recv(recv_buffers[i], 5*recv_counts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Waitall(size,send_req,MPI_STATUSES_IGNORE);

    for(int i=0;i<size;i++){
        if(send_count[i]>0){
            for(int n=0;n<send_count[i];n++){
                int index = n*5;
                PlaneNode *removed = seek_plane(list,(int)send_buffers[i][index]);
                remove_plane(list,removed);
            }
        }
    }

    for(int i=0;i<size;i++){
        if(recv_counts[i]>0){
            for(int n=0;n<recv_counts[i];n++){
                int index = n*5;
                int i_idx = get_index_i(recv_buffers[i][index+1], x_max, N);
                int j_idx = get_index_j(recv_buffers[i][index+2], y_max, M);
                int map_index = get_index(i_idx, j_idx, N, M);
                insert_plane(list,recv_buffers[i][index],map_index, rank, recv_buffers[i][index+1],recv_buffers[i][index+2], recv_buffers[i][index+3], recv_buffers[i][index+4]);
            }
        }
    }

    for (int i=0;i<size;i++) {
        if(send_count[i]>0 && send_buffers[i]!=NULL) {
            free(send_buffers[i]);
        }
        if(recv_counts[i]>0 && recv_buffers[i]!=NULL) {
            free(recv_buffers[i]);
        }
    }
    free(send_buffers);
    free(recv_buffers);
    free(send_count);
    free(recv_counts);
    free(counters);
    free(send_req);
    
}

/// TODO
/// Communicate planes using all to all calls with default data types
void communicate_planes_alltoall(PlaneList* list, int N, int M, double x_max, double y_max, int rank, int size, int* tile_displacements)
{
    int *send_count = calloc(size,sizeof(int));
    
    PlaneNode* current = list->head;
    while (current != NULL) {
        int index_i = get_index_i(current->x, x_max, N);
        int index_j = get_index_j(current->y, y_max, M);
        int plane_rank = get_rank_from_indices(index_i, index_j,  N,  M,tile_displacements, size);
        if(plane_rank != rank){
            send_count[plane_rank]+=1;
        }
        current = current->next;

    }
    for(int i=0;i<size;i++){
        send_count[i] *=5;
    }

    int* recv_counts = malloc(size*sizeof(int));
    MPI_Alltoall(send_count,1,MPI_INT,recv_counts,1,MPI_INT, MPI_COMM_WORLD);

    int* recv_disp = calloc(size,sizeof(int));
    int* send_disp = calloc(size,sizeof(int));

    for(int i=1;i<size;i++){
        recv_disp[i] = recv_disp[i-1] + recv_counts[i-1];
        send_disp[i] = send_disp[i-1] + send_count[i-1];
    }

    int total_recv = recv_disp[size-1] + recv_counts[size-1];
    int total_send = send_disp[size-1] + send_count[size-1];

    double* recv_buffers = malloc(total_recv*sizeof(double));
    double* send_buffers = malloc(sizeof(double)*total_send);
    int* index = calloc(size,sizeof(int));

    current = list->head;
    while (current != NULL) {
        int i = get_index_i(current->x, x_max, N);
        int j = get_index_j(current->y, y_max, M);
        int dest_rank = get_rank_from_indices(i, j, N, M, tile_displacements, size);

        if (dest_rank != rank) {
            int idx = send_disp[dest_rank] + 5*index[dest_rank];
            send_buffers[idx] = (double)current->index_plane;
            send_buffers[idx+1] = current->x;
            send_buffers[idx+2] = current->y;
            send_buffers[idx+3] = current->vx;
            send_buffers[idx+4] = current->vy;
            index[dest_rank]++;
        }
        current = current->next;
    }

    MPI_Alltoallv(send_buffers,send_count,send_disp, MPI_DOUBLE, recv_buffers, recv_counts, recv_disp, MPI_DOUBLE, MPI_COMM_WORLD);

    for(int i=0;i<(total_send/5);i++){
        int idx_plane = i*5;
        PlaneNode* plane = seek_plane(list, (int)send_buffers[idx_plane]);
        remove_plane(list,plane);
    }

    for(int i =0;i<(total_recv/5);i++){
        int pos = i*5;
        int i_idx = get_index_i(recv_buffers[pos+1], x_max, N);
        int j_idx = get_index_j(recv_buffers[pos+2], y_max, M);
        int map_index = get_index(i_idx, j_idx, N, M);
        insert_plane(list,recv_buffers[pos], map_index, rank, recv_buffers[pos+1],recv_buffers[pos+2], recv_buffers[pos+3], recv_buffers[pos+4]);
    }

    free(send_buffers);
    free(recv_buffers);
    free(send_count);
    free(recv_counts);
    free(index);
}

typedef struct {
    int index_plane;
    double x;
    double y;
    double vx;
    double vy;
} MinPlaneToSend;

/// TODO
/// Communicate planes using all to all calls with custom data types
void communicate_planes_struct_mpi(PlaneList* list, int N, int M, double x_max, double y_max, int rank, int size,int* tile_displacements)
{

    MPI_Datatype planetype;
    int blockcounts[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_INT,MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[5];

    offsets[0] = offsetof(MinPlaneToSend,index_plane);
    offsets[1] = offsetof(MinPlaneToSend,x);
    offsets[2] = offsetof(MinPlaneToSend,y);
    offsets[3] = offsetof(MinPlaneToSend,vx);
    offsets[4] = offsetof(MinPlaneToSend,vy);

    MPI_Type_create_struct(5,blockcounts,offsets,types,&planetype);
    MPI_Type_commit(&planetype);

    int* send_count = calloc(size,sizeof(int));
    PlaneNode* current = list->head;
    while(current!=NULL){
        int i = get_index_i(current->x,x_max,N);
        int j = get_index_j(current->y,y_max,M);
        int plane_rank = get_rank_from_indices(i,j,N,M,tile_displacements, size);
        if (plane_rank!=rank){
            send_count[plane_rank]++;
        }
        current = current->next;
    }

    int* recv_count = malloc(size*sizeof(int));
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,MPI_COMM_WORLD);

    int* send_disp = calloc(size,sizeof(int));
    int* recv_disp = calloc(size,sizeof(int));
    for (int i=1;i<size;i++) {
        send_disp[i] = send_disp[i-1] + send_count[i-1];
        recv_disp[i] = recv_disp[i-1] + recv_count[i-1];
    }

    int total_send = send_disp[size-1] + send_count[size-1];
    int total_recv = recv_disp[size-1] + recv_count[size-1];
    
    MinPlaneToSend* send_buffer = malloc(total_send*sizeof(MinPlaneToSend));
    MinPlaneToSend* recv_buffer = malloc(total_recv*sizeof(MinPlaneToSend));
    int* index = calloc(size,sizeof(int));

    current = list->head;
    while(current!=NULL){
        int i = get_index_i(current->x, x_max, N);
        int j = get_index_j(current->y, y_max, M);
        int dest_rank = get_rank_from_indices(i, j, N, M, tile_displacements, size);
        if (dest_rank!=rank) {
            int idx = send_disp[dest_rank] + index[dest_rank]++;
            send_buffer[idx].index_plane = current->index_plane;
            send_buffer[idx].x = current->x;
            send_buffer[idx].y = current->y;
            send_buffer[idx].vx = current->vx;
            send_buffer[idx].vy = current->vy;
        }
        current = current->next;
    }

    MPI_Alltoallv(send_buffer,send_count,send_disp,planetype,recv_buffer, recv_count, recv_disp, planetype, MPI_COMM_WORLD);

    for (int i=0;i<total_send;i++) {
        PlaneNode* plane = seek_plane(list, send_buffer[i].index_plane);
        remove_plane(list,plane);
    }

    for (int i=0;i<total_recv;i++) {
        int i_idx = get_index_i(recv_buffer[i].x, x_max, N);
        int j_idx = get_index_j(recv_buffer[i].y, y_max, M);
        int map_index = get_index(i_idx, j_idx, N, M);
        insert_plane(list,recv_buffer[i].index_plane, map_index, rank,recv_buffer[i].x, recv_buffer[i].y, recv_buffer[i].vx, recv_buffer[i].vy);
    }

    free(send_buffer); 
    free(recv_buffer); 
    free(send_count); 
    free(recv_count); 
    free(send_disp); 
    free(recv_disp);
    free(index);
    MPI_Type_free(&planetype);
}

int main(int argc, char **argv) {
    int debug = 1;                      // 0: no debug, 1: shows all planes information during checking
    int N = 0, M = 0;                   // Grid dimensions
    double x_max = 0.0, y_max = 0.0;    // Total grid size
    int max_steps;                      // Total simulation steps
    char* input_file;                   // Input file name
    int check;                          // 0: no check, 1: check the simulation is correct

    int rank, size;

    /// TODO
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int tile_displacements[size+1];
    int mode = 0;
    if (argc == 5) {
        input_file = argv[1];
        max_steps = atoi(argv[2]);
        if (max_steps <= 0) {
            fprintf(stderr, "max_steps needs to be a positive integer\n");
            return 1;
        }
        mode = atoi(argv[3]);
        if (mode > 2 || mode < 0) {
            fprintf(stderr, "mode needs to be a value between 0 and 2\n");
            return 1;
        }
        check = atoi(argv[4]);
        if (check >= 2 || check < 0) {
            fprintf(stderr, "check needs to be a 0 or 1\n");
            return 1;
        }
    }
    else {
        fprintf(stderr, "Usage: %s <filename> <max_steps> <mode> <check>\n", argv[0]);
        return 1;
    }

    PlaneList owning_planes = {NULL, NULL};
    read_planes_mpi(input_file, &owning_planes, &N, &M, &x_max, &y_max, rank, size, tile_displacements);
    PlaneList owning_planes_t0 = copy_plane_list(&owning_planes);

    //print_planes_par_debug(&owning_planes);

    double time_sim = 0., time_comm = 0, time_total=0;

    double start_time = MPI_Wtime();
    int step = 0;
    for (step = 1; step <= max_steps; step++) {
        double start = MPI_Wtime();
        PlaneNode* current = owning_planes.head;
        while (current != NULL) {
            current->x += current->vx;
            current->y += current->vy;
            current = current->next;
        }
        filter_planes(&owning_planes, x_max, y_max);
        time_sim += MPI_Wtime() - start;

        start = MPI_Wtime();
        if (mode == 0)
            communicate_planes_send(&owning_planes, N, M, x_max, y_max, rank, size, tile_displacements);
        else if (mode == 1)
            communicate_planes_alltoall(&owning_planes, N, M, x_max, y_max, rank, size, tile_displacements);
        else
            communicate_planes_struct_mpi(&owning_planes, N, M, x_max, y_max, rank, size, tile_displacements);
        time_comm += MPI_Wtime() - start;
    }
    time_total = MPI_Wtime() - start_time;

    /// TODO Check computational times
    double total_time_sim = 0;
    double total_time_comm = 0;
    double total_time_total = 0;

    MPI_Reduce(&time_sim, &total_time_sim, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_comm, &total_time_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_total, &total_time_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Flight controller simulation: #input %s mode: %d size: %d\n", input_file, mode, size);
        printf("Time simulation:     %.2fs\n", total_time_sim);
        printf("Time communication:  %.2fs\n", total_time_comm);
        printf("Time total:          %.2fs\n", total_time_total);
    }

    if (check ==1)
        check_planes_mpi(&owning_planes_t0, &owning_planes, N, M, x_max, y_max, max_steps, tile_displacements, size, debug);

    MPI_Finalize();
    return 0;
}
