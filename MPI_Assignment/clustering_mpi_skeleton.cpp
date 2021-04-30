//Yang LI
//Student ID: 20750699
//Email: ylikp@connect.ust.hk

#include "clustering.h"

#include "mpi.h"

#include <cassert>
#include <chrono>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm;
  int num_process; // number of processors
  int my_rank;     // my global rank

  comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &num_process);
  MPI_Comm_rank(comm, &my_rank);

  if (argc != 3) {
    std::cerr << "usage: ./clustering_sequential data_path result_path"
              << std::endl;

    return -1;
  }
  std::string dir(argv[1]);
  std::string result_path(argv[2]);

  int num_graphs;
  int *clustering_results = nullptr;
  int *num_cluster_total = nullptr;

  int *nbr_offs = nullptr, *nbrs = nullptr;
  int *nbr_offs_local = nullptr, *nbrs_local = nullptr;

  GraphMetaInfo *info = nullptr;

  // read graph info from files
  if (my_rank == 0) {
    num_graphs = read_files(dir, info, nbr_offs, nbrs);
  }
  auto start_clock = chrono::high_resolution_clock::now();


  // ADD THE CODE HERE
  // Begin my code part *************

  // Broadcast the graph numbers of each process
  MPI_Bcast(&num_graphs, 1, MPI_INT, 0, comm);
  int local_num_graphs = num_graphs / num_process;

  // Calculate the final result size in process 0
  int length_clustering_results_overall = 0;   // The final result size, to be allocated later
  if (my_rank == 0) 
  {
    for (int i = 0; i < num_graphs; i++)
    {
      length_clustering_results_overall += info->num_vertices;
    }
  }

  // Derived Type for <GraphMetaInfo>, scatter this derived type later
  int array_of_block_length[2] = {1,1};
  MPI_Aint edge_addr, vertice_addr, info_addr;
  MPI_Aint array_of_displacement[2];
  MPI_Datatype array_of_types[2] = {MPI_INT, MPI_INT};
  MPI_Datatype graph_type;
  MPI_Get_address(info, &info_addr);
  MPI_Get_address(&(info->num_edges), &edge_addr);
  MPI_Get_address(&(info->num_vertices), &vertice_addr);
  array_of_displacement[0] = edge_addr - info_addr;
  array_of_displacement[1] = vertice_addr - edge_addr;
  MPI_Type_create_struct(2, array_of_block_length, array_of_displacement, array_of_types, &graph_type); // Create the struct type
  MPI_Type_commit(&graph_type);  // Commit the type before using it

  // Scatter <info> to <local_info> using the derived type <graph_type>
  GraphMetaInfo* local_info;
  local_info = (GraphMetaInfo *)calloc(local_num_graphs, sizeof(GraphMetaInfo));
  if (my_rank == 0)
  {
    MPI_Scatter(info, local_num_graphs, graph_type, local_info, local_num_graphs, graph_type, 0, comm);
  }
  else 
  {
    MPI_Scatter(info, local_num_graphs, graph_type, local_info, local_num_graphs, graph_type, 0, comm);
  }
  
  // Calculate local number of vertices and edges for each process
  int local_space_v(0), local_space_e(0);
  for (int i = 0; i < local_num_graphs; i++)
  {
    local_space_v += local_info[i].num_vertices;
    local_space_e += local_info[i].num_edges;
  }

  // Calculate the arrays of <sendcounts> and <displs> to be used in MPI_Scatterv
  // Because they have different size for each process
  // Need <sendcounts> and <displs> for <nbr_offs>, <nbrs> and <clustering_results>
  int *nbr_offs_sendcounts, *nbr_offs_displs, *nbrs_sendcounts, *nbrs_displs, *clustering_results_recvcounts, *clustering_results_displs;
  nbr_offs_sendcounts = (int*)calloc(num_process, sizeof(int));
  nbrs_sendcounts = (int*)calloc(num_process, sizeof(int));
  nbr_offs_displs = (int*)calloc(num_process, sizeof(int));
  nbrs_displs = (int*)calloc(num_process, sizeof(int));
  clustering_results_recvcounts = (int*)calloc(num_process, sizeof(int));
  clustering_results_displs = (int*)calloc(num_process, sizeof(int));

  if (my_rank == 0) // Process 0 has the whole data, so calculate the size in process 0
  {
    int num_v(0), num_e(0), dis_v(0), dis_e(0);
    int j = 0; // The process number
    for (int i = 0; i < num_graphs; i++)  
    {
      num_v += info[i].num_vertices;
      num_e += info[i].num_edges;
      if ( (i+1) % local_num_graphs == 0 && i != 0) // Every <local_num_graphs>, store the counts and displs
      {
        nbr_offs_sendcounts[j] = num_v + local_num_graphs;  // <nbrs> and <nbr_offs> need 1 more int space for each graph
        nbrs_sendcounts[j] = num_e + local_num_graphs;
        clustering_results_recvcounts[j] = num_v;

        dis_v += num_v;
        dis_e += num_e;
        if (j + 1 < num_process)
        {
          nbr_offs_displs[j+1] = dis_v + local_num_graphs * (j+1);
          nbrs_displs[j+1] = dis_e + local_num_graphs * (j+1);
          clustering_results_displs[j+1] = dis_v;
        }
        num_v = 0;
        num_e = 0;
        j++;
      }
    }
  }

  // Broadcast the sendcounts and displs from process 0 
  MPI_Bcast(nbr_offs_sendcounts, num_process, MPI_INT, 0, comm);
  MPI_Bcast(nbr_offs_displs, num_process, MPI_INT, 0, comm);
  MPI_Bcast(nbrs_sendcounts, num_process, MPI_INT, 0, comm);
  MPI_Bcast(nbrs_displs, num_process, MPI_INT, 0, comm);
  MPI_Bcast(clustering_results_recvcounts, num_process, MPI_INT, 0, comm);
  MPI_Bcast(clustering_results_displs, num_process, MPI_INT, 0, comm);

  // Allocate space for local nbr_offs, nbrs with the previous result 
  // And then use MPI_Scatterv to distribute the values
  nbr_offs_local = (int*)calloc(local_space_v + local_num_graphs, sizeof(int));
  nbrs_local = (int*)calloc(local_space_e + local_num_graphs, sizeof(int));
  if (my_rank == 0)
  {
    MPI_Scatterv(nbr_offs, nbr_offs_sendcounts, nbr_offs_displs, MPI_INT, nbr_offs_local, local_space_v + 3, MPI_INT, 0, comm);
    MPI_Scatterv(nbrs, nbrs_sendcounts, nbrs_displs, MPI_INT, nbrs_local, local_space_e + 3, MPI_INT, 0, comm);
  }
  else
  {
    MPI_Scatterv(nbr_offs, nbr_offs_sendcounts, nbr_offs_displs, MPI_INT, nbr_offs_local, local_space_v + 3, MPI_INT, 0, comm);
    MPI_Scatterv(nbrs, nbrs_sendcounts, nbrs_displs, MPI_INT, nbrs_local, local_space_e + 3, MPI_INT, 0, comm);
  }

  // Allocate space for local results
  int* local_clustering_results;
  int* local_num_cluster;
  local_clustering_results = (int *)calloc(local_space_v, sizeof(int));
  local_num_cluster = (int*)calloc(local_num_graphs, sizeof(int));

  // Do graph clustering locally
  int* cus_clustering_results = local_clustering_results;
  for (int i = 0; i < local_num_graphs; i++)
  {
    GraphMetaInfo cus_info = local_info[i]; 
    int cus_num_cluster = clustering(cus_info, nbr_offs_local, nbrs_local, cus_clustering_results);
    local_num_cluster[i] = cus_num_cluster;
    //printf("num cluster in graph %d : %d\n", i + my_rank * 3, cus_num_cluster);
      
    cus_clustering_results += (cus_info.num_vertices);
    nbr_offs_local += (cus_info.num_vertices + 1);
    nbrs_local += (cus_info.num_edges + 1);       
  }

  // Use MPI_Gatherv to collect the local results to process 0
  if (my_rank == 0) {
    num_cluster_total = (int*)calloc(num_graphs, sizeof(int));
    clustering_results = (int*)calloc(length_clustering_results_overall, sizeof(int));
    MPI_Gather(local_num_cluster, local_num_graphs, MPI_INT, num_cluster_total, local_num_graphs, MPI_INT, 0, comm);
    MPI_Gatherv(local_clustering_results, local_space_v, MPI_INT, clustering_results, clustering_results_recvcounts, clustering_results_displs,  MPI_INT, 0, comm);
  }
  else {
    MPI_Gather(local_num_cluster, local_num_graphs, MPI_INT, num_cluster_total, local_num_graphs, MPI_INT, 0, comm);
    MPI_Gatherv(local_clustering_results, local_space_v, MPI_INT, clustering_results, clustering_results_recvcounts, clustering_results_displs, MPI_INT, 0, comm);
  }
  // Finalize my code part *************

  MPI_Barrier(comm);
  auto end_clock = chrono::high_resolution_clock::now();

  // 1) print results to screen
  if (my_rank == 0) {
    for (size_t i = 0; i < num_graphs; i++) {
      printf("num cluster in graph %d : %d\n", i, num_cluster_total[i]);
    }
    fprintf(stderr, "Elapsed Time: %.9lf ms\n",
            chrono::duration_cast<chrono::nanoseconds>(end_clock - start_clock)
                    .count() /
                pow(10, 6));
  }

  // 2) write results to file
  if (my_rank == 0) {
    int *result_graph = clustering_results;
    for (int i = 0; i < num_graphs; i++) {
      GraphMetaInfo info_local = info[i];
      write_result_to_file(info_local, i, num_cluster_total[i], result_graph,
                           result_path);

      result_graph += info_local.num_vertices;
    }
  }

  MPI_Type_free(&graph_type);
  MPI_Finalize();

  if (my_rank == 0) {
    free(num_cluster_total);
  }

  return 0;
}
