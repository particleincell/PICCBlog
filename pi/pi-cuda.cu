#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

using namespace std;

constexpr int BLOCK_DIM = 16;

// determines if this node (pixel) is inside the circle
// result is stored in a [16*16] array
// thread 0 then computes the number of "in" nodes (value from 0 to 16*16)
__global__ void flagKernel(unsigned *block_counts) {
	bool __shared__ ins[BLOCK_DIM*BLOCK_DIM];
	// compute our coordinate in the global grid
	unsigned i = blockIdx.x*blockDim.x + threadIdx.x; // my i
	unsigned j = blockIdx.y*blockDim.y + threadIdx.y; // my j
	unsigned Ni = gridDim.x*blockDim.x;   // total number of nodes in x
	unsigned Nj = gridDim.y*blockDim.y;   // total number of nodex in y

	//get 1D index from i,j, u=j*ni+i
	unsigned u = threadIdx.y*blockDim.x + threadIdx.x;

	float x = i/(float)Ni;     // compute x in [0,1)
	float y = j/(float)Nj;     // y in [0,1)
	if (x*x+y*y<=1) ins[u] = true;  // check if in the circle
	else ins[u] = false;

	// wait for all threads in the block to finish
	__syncthreads();

	// let the first thread in the block add up "ins"
	if (u==0) {
		unsigned count = 0;
		for (int i=0;i<blockDim.x*blockDim.y;i++)
		  if (ins[u]) count++;

		// flattened index for the block, u=j*ni+i
		int block_u = blockIdx.y*gridDim.x+blockIdx.x;

		// store the sum in global memory
		block_counts[block_u] = count;
	}
}

// this kernel adds up block-level sums to the global sum
// this could be further optimized by splitting up the sum over threads
__global__ void addKernel(dim3 numBlocks, unsigned *block_counts, unsigned long *glob_count) {
	// compute total number of blocks
	unsigned N = numBlocks.x*numBlocks.y;
	unsigned long sum = 0;
	for (int i=0;i<N;i++)
		sum+=block_counts[i];

	// store result in global memory
	*glob_count = sum;
}


int main() {
  // grab starting time
  auto time_start = chrono::high_resolution_clock::now();

  // figure out how many samples I should process
  size_t N = BLOCK_DIM*1000;    // grid size

  // figure out our grid size
  dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM);
  dim3 numBlocks(N / threadsPerBlock.x, N / threadsPerBlock.y);

  // allocate memory on the GPU
  unsigned *block_counts;
  cudaMalloc((void**)&block_counts, numBlocks.x*numBlocks.y*sizeof(unsigned));

  unsigned long *N_in_gpu;  // GPU variable to hold the total N_in
  unsigned long N_in;		// CPU variable to hold this data
  cudaMalloc((void**)&N_in_gpu, sizeof(N_in));

  // launch the kernel to flag nodes, each block has BLOCK_DIM*BLOCK_DIM threads
  flagKernel<<<numBlocks, threadsPerBlock>>>(block_counts);

  // launch kernel to add up per-block "in" counts
  addKernel<<<1, 1>>>(numBlocks, block_counts, N_in_gpu);

  // transfer N_in from the GPU to the CPU
  cudaMemcpy(&N_in, N_in_gpu, sizeof(N_in), cudaMemcpyDeviceToHost);

  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;

  // compute pi and show the result on rank 0 (root) using the global data
  size_t N_tot = N*N;
  double pi = 4*N_in/(double)N_tot;
  cout<<"Using a "<<N<<"x"<<N<<" grid ("<<N_tot<<" samples), pi is "<<pi
	  <<" in "<<setprecision(3)<<time_delta.count()<<" seconds"<<endl;

  // be a good neighbor and free memory
  cudaFree(block_counts);
  cudaFree(N_in_gpu);

  return 0;
}

