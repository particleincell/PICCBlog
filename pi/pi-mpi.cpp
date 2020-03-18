#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <sstream>
#include <mpi.h>

using namespace std;

// object for sampling random numbers
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);} 

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

Rnd rnd;		// instantiate array of Rnd objects

int main(int num_args, char**args) {
  MPI_Init(&num_args, &args);      // initialize MPI

  // grab starting time
  auto time_start = chrono::high_resolution_clock::now();

  // figure out what rank I am and the total number of processes
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

  // write to a stringstream first to prevent garbled output
  /*stringstream ss;
    ss<<"I am "<<mpi_rank<<" of "<<mpi_size<<endl;
    cout<<ss.str();*/

  // figure out how many samples I should process
//  unsigned N_goal = 100000000;    
  unsigned N_goal = (16*1000)*(16*1000);    

  unsigned N_tot =  N_goal / mpi_size;
  
  // correct on last rank to get the correct total
  if (mpi_rank == mpi_size-1) N_tot = N_goal - (N_tot*(mpi_size-1));
	 
  // perform the computation
  unsigned N_in = 0;	
  for (size_t s=0;s<N_tot;s++) {
     double x = rnd();    
 	 double y = rnd();    
	 if (x*x+y*y<=1) N_in++;
  }   

  // sum up N_in and N_tot across all ranks and send the sum to rank 0
  unsigned N_tot_glob;
  unsigned N_in_glob;
  MPI_Reduce(&N_in,&N_in_glob,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&N_tot,&N_tot_glob,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);

  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  
  // compute pi and show the result on rank 0 (root) using the global data
  if (mpi_rank==0) {
    double pi = 4*N_in_glob/(double)N_tot_glob; 
    stringstream ss;
    ss<<"Using "<<N_tot_glob<<" samples and "<<mpi_size<<" processes, pi is "<<pi
	  <<" in "<<setprecision(3)<<time_delta.count()<<" seconds"<<endl;
    cout<<ss.str();
  }
  
  MPI_Finalize();     // indicate clean exit
  return 0;
}
