#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <vector>

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

vector<Rnd> rnds;		// instantiate array of Rnd objects

// function to launch in parallel, we could alternatively define a class
// result is stored in N_in[thread_id]
void Worker(int thread_id, size_t N_total, size_t *N_ins) {
    // set references to thread-specific items
	size_t &N_in = N_ins[thread_id];	// result
	Rnd &rnd = rnds[thread_id];			// our dedicated generator
	N_in = 0;	
  	for (size_t s=0;s<N_total;s++) {
    	double x = rnd();    
 		double y = rnd();    
		if (x*x+y*y<=1) N_in++;
  	}   
  }

int main(int num_args, char**args) {
  // get maximum number of supported concurrent threads
  int num_threads = thread::hardware_concurrency();
  if (num_args>1) num_threads = atoi(args[1]); // get from command line

  // generate an array of random number generators
  rnds.resize(num_threads);

  // grab starting time
  auto time_start = chrono::high_resolution_clock::now();

  size_t N_goal = 100000000;
  size_t N_total = 0;

  // allocate memory for results
  size_t *N_in_arr = new size_t[num_threads];

  // container to store thread objects
  vector<thread> threads;

  // create threads
  for (int i=0;i<num_threads;i++) {
     size_t N_per_thread =  N_goal / num_threads;
	 // correct on last one to get the correct total
     if (i==num_threads-1) N_per_thread = N_goal - (N_per_thread*(num_threads-1));
     N_total += N_per_thread;
	 
     // launch Worker(i,N_per_thread,N_in_arr) in parallel
     threads.emplace_back(Worker,i,N_per_thread,N_in_arr);
  }
  
  // wait for workers to finish
  for (thread &th:threads) 
   th.join();
  
  // sum up thread computed results
  size_t N_in = 0;
  for (int i=0;i<num_threads;i++) {     
     N_in += N_in_arr[i];
  }

  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  
  double pi = 4*N_in/(double)N_total;  // cast to double to perform floating point division
  cout<<"Using "<<N_total<<" samples and "<<num_threads<<" threads, pi is "<<pi
	  <<" in "<<setprecision(3)<<time_delta.count()<<" seconds"<<endl;

  delete[] N_in_arr;
  
  return 0;
}
