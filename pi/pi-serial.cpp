#include <iostream>
#include <iomanip>
#include <random>
#include <math.h>

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

Rnd rnd;		// instantiate the Rnd object

int main() {
  size_t N_total = 10000000;
  size_t N_in = 0;   

  for (size_t s=0;s<N_total;s++) {
    double x = rnd();       // pick random x in [0,1)
 	double y = rnd();       // pick random y in [0,1)
	if (x*x+y*y<=1) N_in++;	// increment counter if inside the circle
  }

  double pi = 4*N_in/(double)N_total;  // cast to double to perform floating point division
  double error = 100*abs(pi/acos(-1.0)-1); // (our-real)/real
  cout<<"Using "<<N_total<<" samples, pi is "<<pi<<" ("<<setprecision(2)<<error<<"% error)"<<endl;
  
  return 0;
}
