/*
 * Name: LI Yang
 * Student id: 20750699
 * ITSC email: ylikp@connect.ust.hk
 *
 * Please only change this file and do not change any other files.
 * Feel free to change/add any helper functions.
 *
 * (problem)COMPILE: g++ -lstdc++ -std=c++11 -lpthread clustering_pthread_skeleton.cpp -main.cpp -o pthread
 * (ok)COMPILE: g++ -lstdc++ -std=c++11 clustering_pthread_skeleton.cpp main.cpp -o pthread -lpthread 
 * RUN:     ./pthread <path> <epsilon> <mu> <num_threads>
 */

#include <pthread.h>
#include "clustering.h"
 
long long num_vs_G;
long long num_es_G;
bool* pivots_G = nullptr;
int* num_sim_nbrs_G = nullptr;
int** sim_nbrs_G = nullptr;
int* nbr_offs_G = nullptr;
int* nbrs_G = nullptr;
bool* visited_G = nullptr;
int* clustering_result = nullptr;
float epsilon_G;
int mu_G;
pthread_mutex_t mutex;
pthread_cond_t cond_var;
pthread_rwlock_t rwlock;
int counter = 0;

struct AllThings{
    int num_threads;
    int my_rank;

    AllThings(int inum_threads, int imy_rank){
        num_threads = inum_threads;
        my_rank = imy_rank;
    };
};

void expansion(int cur_id, int num_clusters, int *num_sim_nbrs, int **sim_nbrs,
               bool *visited, bool *pivots, int *cluster_result) 
{
  for (int i = 0; i < num_sim_nbrs[cur_id]; i++) 
  {
    int nbr_id = sim_nbrs[cur_id][i];

    if ((pivots[nbr_id])&&(!visited[nbr_id]))
    {
        visited[nbr_id] = true;
        //printf("-------------------------------------ready to change result %d to be %d \n", nbr_id, num_clusters);
        //fflush(stdout);
        pthread_rwlock_wrlock(&rwlock);
        if ( cluster_result[nbr_id] > num_clusters || cluster_result[nbr_id] == -1 )
        {
            //printf("------------------------------------------------------------------------change result %d to be %d \n", nbr_id, num_clusters);
            //fflush(stdout);
            cluster_result[nbr_id] = num_clusters;  
        }
        pthread_rwlock_unlock(&rwlock);

        expansion(nbr_id, num_clusters, num_sim_nbrs, sim_nbrs, visited, pivots,
                cluster_result);
    }
  }
}


void *parallel(void* allthings){
    AllThings *all = (AllThings *) allthings;

    printf("Hello from %d of %d\n", all->my_rank, all->num_threads);
    long my_rank = (long)(all->my_rank);
    int thread_count = (long)(all->num_threads);
    long long local_num_vs = ( (num_vs_G - 1) / thread_count ) + 1;

    int my_first = my_rank * local_num_vs;
    int last = (my_rank + 1) * local_num_vs - 1;
    int my_last = (last < num_vs_G ? last : (num_vs_G - 1));
    //printf("last is %d, num_vs_G is %lld \n", last, num_vs_G);
    printf("my rank is %ld, my_first is %d, my_last is %d --- \n", my_rank, my_first, my_last);
    fflush(stdout);
    for (int i = my_first; i < my_last + 1; i++) 
    {
        int *left_start = &nbrs_G[nbr_offs_G[i]];
        int *left_end = &nbrs_G[nbr_offs_G[i + 1]];
        int left_size = left_end - left_start;

        sim_nbrs_G[i] = new int[left_size];
    // loop over all neighbors of i
        for (int *j = left_start; j < left_end; j++) 
        {
            int nbr_id = *j;

            int *right_start = &nbrs_G[nbr_offs_G[nbr_id]];
            int *right_end = &nbrs_G[nbr_offs_G[nbr_id + 1]];
            int right_size = right_end - right_start;

            // compute the similarity
            int num_com_nbrs = get_num_com_nbrs(left_start, left_end, right_start, right_end);

            float sim = (num_com_nbrs + 2) / std::sqrt((left_size + 1.0) * (right_size + 1.0));

            if (sim > epsilon_G) {
                sim_nbrs_G[i][num_sim_nbrs_G[i]] = nbr_id;
                num_sim_nbrs_G[i]++;
            }
        }

        if (num_sim_nbrs_G[i] > mu_G) pivots_G[i] = true;
    }

    pthread_mutex_lock(&mutex);
    counter ++;
    if (counter == thread_count)
    {
        counter = 0;
        pthread_cond_broadcast(&cond_var);
    }
    else
    {
        while (pthread_cond_wait(&cond_var, &mutex) != 0);
    }
    pthread_mutex_unlock(&mutex);
    //stage 2

    bool* local_visited = new bool[num_vs_G]();

    for (int i = my_first; i < my_last +1; i++) 
    {
        printf("from thread %ld, pivots %d is %d \n", my_rank, i, pivots_G[i]);
        fflush(stdout);

        if (!pivots_G[i] || local_visited[i]) continue;

        visited_G[i] = true;
        local_visited[i] = true;

       
        if ( clustering_result[i] > i || clustering_result[i] == -1 )
        {
            //pthread_rwlock_wrlock(&rwlock);
            //printf("-------------------------------------external change %d, result %d to be %d \n", clustering_result[i], i, i);
            //fflush(stdout);
            clustering_result[i] = i;
            expansion(i, i, num_sim_nbrs_G, sim_nbrs_G, local_visited, pivots_G, clustering_result);
        }
        else continue;
    }
  
    return 0;
}

int *scan(float epsilon, int mu, int num_threads, int num_vs, int num_es, int *nbr_offs, int *nbrs){
    long thread;
    pthread_t* thread_handles = (pthread_t*) malloc(num_threads*sizeof(pthread_t));
    //int *cluster_result = new int[num_vs];

    num_vs_G = num_vs;
    num_es_G = num_es;
    pivots_G = new bool[num_vs]();
    num_sim_nbrs_G = new int[num_vs]();
    sim_nbrs_G = new int*[num_vs];
    nbr_offs_G = nbr_offs;
    nbrs_G = nbrs;
    epsilon_G = epsilon;
    mu_G = mu;

    visited_G = new bool[num_vs]();
    clustering_result = new int[num_vs];
    std::fill(clustering_result, clustering_result + num_vs_G, -1);

    pthread_mutex_init(&mutex, NULL);
    pthread_rwlock_init(&rwlock, NULL);

    for (thread=0; thread < num_threads; thread++)
        pthread_create(&thread_handles[thread], NULL, parallel, (void *) new AllThings(
                num_threads, thread));
    for (thread=0; thread < num_threads; thread++)
        pthread_join(thread_handles[thread], NULL);

    pthread_mutex_destroy(&mutex);
    pthread_rwlock_destroy(&rwlock);


/*
    for (int i = 0; i < num_vs_G; i++)
    {
        printf("%d th of pivots is %d \n", i, pivots_G[i]);
        fflush(stdout);
    }
*/
    delete[] pivots_G;
    delete[] num_sim_nbrs_G;
    for (auto i = 0; i < num_vs; i++)
    {
        delete[] sim_nbrs_G[i];
    }
    delete[] sim_nbrs_G;
    return clustering_result;
}



