#include <stdio.h>
#include <stdlib.h>

__device__ int d_value;
__global__ void test_Kernel()
{
	int threadID = threadIdx.x;
	d_value = 1;
	printf("threadID %-3d d_value%3d\n",threadID,d_value);
}
int main()
{
	int h_value = 0;
	test_Kernel<<<1,2>>>();
	cudaMemcpyFromSymbol(&h_value,d_value,
		sizeof(int),0,cudaMemcpyDeviceToHost);
	
	printf("Output from host: %d\n",h_value);
	return 0;
}

