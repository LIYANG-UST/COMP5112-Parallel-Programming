/* 
 * LI Yang
 * 20750699
 * ylikp@connect.ust.hk
 *
 * COMPILE: nvcc -std=c++11 clustering_cuda_skeleton.cu clustering_impl.cpp main.cpp -o cuda
 * RUN:     ./cuda <path> <epsilon> <mu> <num_blocks_per_grid> <num_threads_per_block>
 */

#include <iostream>
#include "clustering.h"

// Define variables or functions here
__device__ int get_num_com_nbrs(int *nbrs, int left_start, int left_end, int right_start, int right_end) 
{
    int left_pos = left_start, right_pos = right_start, num_com_nbrs = 0;

    while (left_pos < left_end && right_pos < right_end) {
        if (nbrs[left_pos] == nbrs[right_pos]) {
            num_com_nbrs++;
            left_pos++;
            right_pos++;
        } else if (nbrs[left_pos] < nbrs[right_pos]) {
            left_pos++;
        } else {
            right_pos++;
        }
    }
    return num_com_nbrs;
}

void expansion(int cur_id, int num_clusters, int *num_sim_nbrs, int **sim_nbrs,
    bool *visited, bool *pivots, int *cluster_result)
{
    for (int i = 0; i < num_sim_nbrs[cur_id]; i++) 
    {
        int nbr_id = sim_nbrs[cur_id][i];
        if ((pivots[nbr_id])&&(!visited[nbr_id]))
        {
            visited[nbr_id] = true;
            cluster_result[nbr_id] = num_clusters;
            expansion(nbr_id, num_clusters, num_sim_nbrs, sim_nbrs, visited, pivots,
                cluster_result);
        }
    }
}


__global__ void stage_1(int num_vs_G, int* Device_nbr_offs, int* Device_nbrs, 
                        int num_blocks_per_grid, int num_threads_per_block, 
                        bool* pivots, int* num_sim_nbrs, int* sim_nbrs,
                        int* size_index, float epsilon, int mu)
{
  
    int my_thread_rank = blockDim.x * blockIdx.x + threadIdx.x;
    int num_threads = num_blocks_per_grid * num_threads_per_block;
    /*
    //Change to use the coalesced access (faster than the way used in ass2)
    int local_num_vs = ( (num_vs_G - 1) / num_threads ) + 1;
    int my_first = my_thread_rank * local_num_vs;
    int last = (my_thread_rank + 1) * local_num_vs - 1;
    int my_last = (last < num_vs_G ? last : (num_vs_G - 1));
    */
    //printf("my_rank is %d, my_first is %d, my_last is %d \n", my_thread_rank, my_first, my_last);

    for (int i = my_thread_rank; i < num_vs_G; i+= num_threads) 
    {
        int left_start = Device_nbr_offs[i];
        int left_end = Device_nbr_offs[i + 1];
        int left_size = left_end - left_start;
        //printf("device:  i %d, left_start is %d, left_end is %d, left_size is %d \n", i, left_start, left_end, left_size);
        //sim_nbrs[i] = new int[left_size];
        // loop over all neighbors of i
        for (int j = left_start; j < left_end; j++) {
            int nbr_id = Device_nbrs[j];

            int right_start = Device_nbr_offs[nbr_id];
            int right_end = Device_nbr_offs[nbr_id + 1];
            int right_size = right_end - right_start;

            // compute the similarity
            int num_com_nbrs = get_num_com_nbrs(Device_nbrs, left_start, left_end, right_start, right_end);

            float sim = (num_com_nbrs + 2) / std::sqrt((left_size + 1.0) * (right_size + 1.0));

            if (sim > epsilon) {
                //sim_nbrs[i][num_sim_nbrs[i]] = nbr_id;
                sim_nbrs[size_index[i]+num_sim_nbrs[i]] = nbr_id;
                num_sim_nbrs[i]++;
            }
        }
        if (num_sim_nbrs[i] > mu) pivots[i] = true;
    }
    __syncthreads();
}


/* 
// how to implement stage 2 on GPU without the read/write lock ??

__global__ void stage_2(int num_vs_G, int num_blocks_per_grid, int num_threads_per_block, 
                        bool* pivots, int* num_sim_nbrs, int* sim_nbrs,
                        bool* visited, int* cluster_result, int* num_clusters,
                        int* size_index)
{
    int my_thread_rank = blockDim.x * blockIdx.x + threadIdx.x;
    int num_threads = num_blocks_per_grid * num_threads_per_block;

    for (int i = my_thread_rank; i < num_vs_G; i += num_threads) 
    {
        if (!pivots[i] || visited[i]) continue;

        if ( cluster_result[i] > i || cluster_result[nbr_id] == -1 )
        {
            visited[i] = true;
            cluster_result[i] = i;
            expansion(i, i, num_sim_nbrs, sim_nbrs, visited, pivots, cluster_result, size_index);
            (*num_clusters)++;
        }
        else continue;

    }
    __syncthreads();
}
*/



// Host code
void cuda_scan(int num_vs, int num_es, int *nbr_offs, int *nbrs,
        float epsilon, int mu, int num_blocks_per_grid, int num_threads_per_block,
        int &num_clusters, int *cluster_result) 
    {
        // pivots for host and device
        bool* Device_pivots, *Host_pivots;
        cudaMalloc((void**)&Device_pivots, sizeof(int) * num_vs );
        Host_pivots = new bool[num_vs]();
        // sim_nbrs for host and device. 2D array to 1D array, use 1D array on GPU
        int* Host_num_sim_nbrs;
        Host_num_sim_nbrs = new int[num_vs]();
        int** Host_sim_nbrs = new int*[num_vs]();
        int* Host_num_nbrs = new int[num_vs]();
        int* Host_size_index = new int[num_vs]();
        int left_start, left_end, left_size, size_index(0);
        for (int i = 0; i < num_vs; i++)
        {
            left_start = nbr_offs[i];
            left_end = nbr_offs[i+1];
            left_size = left_end - left_start;
            Host_sim_nbrs[i] = new int[left_size];
            Host_num_nbrs[i] = left_size;
            Host_size_index[i] = size_index;
            size_index += left_size;
        }

        int* Device_num_sim_nbrs, *Device_sim_nbrs;
        cudaMalloc((void**)&Device_num_sim_nbrs, sizeof(int) * num_vs );
        cudaMalloc((void**)&Device_sim_nbrs, sizeof(int) * size_index );
        // Pass the size index to GPU
        int* Device_size_index;
        cudaMalloc((void**)&Device_size_index, sizeof(int) * num_vs );
        cudaMemcpy(Device_size_index, Host_size_index, sizeof(int) * num_vs, cudaMemcpyHostToDevice);
        // Malloc and copy nbr_offs and nbrs to device. <global memory>
        int *Device_nbr_offs, *Device_nbrs;
        cudaMalloc((void**)&Device_nbr_offs, sizeof(int) * (num_vs + 1));
        cudaMalloc((void**)&Device_nbrs, sizeof(int) * (num_es + 1));
	    cudaMemcpy(Device_nbr_offs, nbr_offs, sizeof(int) * (num_vs + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(Device_nbrs, nbrs, sizeof(int) * (num_es + 1), cudaMemcpyHostToDevice);
 

        
        stage_1<<<num_blocks_per_grid, num_threads_per_block>>>(num_vs, Device_nbr_offs, Device_nbrs, 
                                                                num_blocks_per_grid, num_threads_per_block, 
                                                                Device_pivots, Device_num_sim_nbrs, Device_sim_nbrs, 
                                                                Device_size_index, epsilon, mu);
        cudaDeviceSynchronize();
        // Pass back the pivots results and the sim_nbrs results for stage 2                    
        cudaMemcpy(Host_pivots, Device_pivots, sizeof(bool) * num_vs, cudaMemcpyDeviceToHost);
        cudaMemcpy(Host_num_sim_nbrs, Device_num_sim_nbrs, sizeof(int) * num_vs, cudaMemcpyDeviceToHost);
        for (int i = 0; i < num_vs; i++)
        {
            cudaMemcpy(Host_sim_nbrs[i], Device_sim_nbrs + Host_size_index[i], sizeof(int) * Host_num_nbrs[i], cudaMemcpyDeviceToHost);
        }
        
        // Stage 2
        bool* visited = new bool[num_vs]();
        for (int i = 0; i < num_vs; i++)
        {
            if (!Host_pivots[i] || visited[i]) continue;
      
            visited[i] = true;
            cluster_result[i] = i;
            expansion(i, i, Host_num_sim_nbrs, Host_sim_nbrs, visited, Host_pivots, cluster_result);

            num_clusters++;
        }
        /*
        bool* Device_visited;
        cudaMalloc((void**)&Device_visited, sizeof(int) * num_vs);
        int* Device_cluster_result;
        cudaMalloc((void**)&Device_cluster_result, sizeof(int) * num_vs);
        int* Device_num_clusters;
        cudaMalloc((void**)&Device_num_clusters, sizeof(int));
        cudaMemcpy(Device_cluster_result, cluster_result, sizeof(int) * num_vs);
        stage_2<<<num_blocks_per_grid, num_threads_per_block>>>(num_vs, num_blocks_per_grid, num_threads_per_block, 
                                                                Device_pivots, Device_num_sim_nbrs, Device_sim_nbrs,
                                                                Device_visited, Device_cluster_result, Device_num_clusters,
                                                                Device_size_index);
        cudaMemcpy(cluster_result, Device_cluster_result, sizeof(int) * num_vs, cudaMemcpyDeviceToHost);
        cudaMemcpy(&num_clusters, Device_num_clusters, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(Host_pivots, Device_pivots, sizeof(bool) * num_vs, cudaMemcpyDeviceToHost);
        */

        for(int i = 0; i < num_vs; i++)
        {
            delete[] Host_sim_nbrs[i];
        }
        delete[] Host_sim_nbrs;
        delete[] Host_num_nbrs;
        delete[] Host_pivots;
        delete[] Host_num_sim_nbrs;
        delete[] visited;
        delete[] Host_size_index;
        cudaFree(Device_nbr_offs);
        cudaFree(Device_nbrs);
        cudaFree(Device_num_sim_nbrs);
        cudaFree(Device_sim_nbrs);
        cudaFree(Device_pivots);
        cudaFree(Device_size_index);
        
        // Fill in the cuda_scan function here
}
