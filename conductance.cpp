/*
Simulation of molecular conductance through a cylindrical tube
Tube diameter is 1. The simulation runs for 1,000,000 time steps.

syntax: ./conductance L num_threads np
	L			tube length (1)
	num_threads number of threads to use (default to max available)
	np			approximate number of particles to generate (1e7)

More info at https://www.particleincell.com/2016/conductance/
*/
#include <iostream>
#include <iomanip>
#include <thread>
#include <random>
#include <list>

using namespace std;

const double PI = acos(-1.0);

/* vector helper functions */
namespace vec
{
	// r = a
	inline void copy(double r[3], double a[3]) {r[0]=a[0];r[1]=a[1];r[2]=a[2];}

	// r = a + b*s
	inline void mult_and_add(double r[3], double a[3], double b[3], double s) {r[0]=a[0]+b[0]*s;r[1]=a[1]+b[1]*s;r[2]=a[2]+b[2]*s;}

	// r = a *s;
	inline void mult(double r[3], double a[3], double s) {r[0]=a[0]*s;r[1]=a[1]*s;r[2]=a[2]*s;}

	// r = a cross b
	inline void cross(double r[3], double a[3], double b[3]) 
	{
		r[0]=a[1]*b[2]-a[2]*b[1]; 
		r[1]=a[2]*b[0]-a[0]*b[2]; 
		r[2]=a[0]*b[1]-a[1]*b[0];
	}

	// r = {0,0,0}
	inline void zero(double r[3]) {r[0]=r[1]=r[2]=0;}

	// returns vector magnitude
	inline double mag(double a[3]) {return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}
	
	// r = a/mag(a)*/
	inline void unit(double r[3], double a[3]) {double vec_mag=mag(a); if (vec_mag>0) mult(r,a,1.0/mag(a)); else zero(r);}
}

/***************************************************************************/
struct Particle
{
    double pos[3];
    double vel[3];   

	/*constructor*/
    Particle (double pos[3], double vel[3])
    {
        for (int i=0;i<3;i++) this->pos[i]=pos[i];
        for (int i=0;i<3;i++) this->vel[i]=vel[i];
    }
};


/*conductance class*/
class ConductanceXYZ 
{
public:
	const double R=0.5;      //cylinder radius
	const double L;			 //cylinder length
	
	const long np_per_ts;	//number of particles to create per time step
	const int num_ts;		//number of timesteps to run

	const double vdrift = 5;	//injection velocity
	list<Particle*> particles;	//list of particles

	const double dt=1e-2;	//simulation time step size
	int ts;		//current time step
	
	int inlet = 0;	/*number of particles crossing inlet*/
	int outlet = 0; /*number of particles crossing outlet*/
	
	bool finished = false;	//set true once loop finishes

	/*constructor, initializes parameters and seeds random number generator*/
	ConductanceXYZ(double L, long np_per_ts, int num_ts):
		gen((unsigned int)chrono::system_clock::now().time_since_epoch().count()),
		L(L), np_per_ts(np_per_ts), num_ts(num_ts)	{
		thr = new thread(start, this);
	}

	thread *thr;
	/*thread start*/
	static void start(ConductanceXYZ *p)   {p->mainLoop();}

	/*random number generator*/
	mt19937 gen;    	
	double rnd() {
		uniform_real_distribution<double> rnd_dist(0.0,1.0);
		return rnd_dist(gen);
	}

	/*main loop, runs for num_ts time steps*/
	void mainLoop()
	{
		for (ts=0;ts<num_ts;ts++)
		{
			InjectParticles();
			MoveParticles();
		}

		finished = true;
	}

	/*injects np_per_ts particles at z=0 with cosine distribution about z axis*/
	void InjectParticles()
	{
		double source_n[3] = {0,0,1};	//source normal direction
		double tang1[3] = {1,0,0};		//tangent 1
		double tang2[3] = {0,1,0};		//tangent 2

		for (long p=0;p<np_per_ts;p++)
		{
			double pos[3];
			double vel[3];
			double rr;
        
			//select random position between [-R:R]x[-R:R], accept if in tube
			do {
				pos[0]=R*(-1+2*rnd());
				pos[1]=R*(-1+2*rnd());
          
				rr= pos[0]*pos[0]+pos[1]*pos[1];
				if (rr<=R*R) break;
			} while(true);
        
			pos[2]=0;	//z=0
        
			//velocity
			double n[3];		
			lambertianVector(n,source_n,tang1,tang2);

			//multiply by drift velocity
			vel[0]= n[0]*vdrift;
			vel[1]= n[1]*vdrift;
			vel[2]= n[2]*vdrift;
		
			//add to dynamic storage
			particles.push_back(new Particle(pos,vel));
		}    
	}

	/*moves particles, reflecting ones hitting the tube wall*/
	void MoveParticles()
	{
	   double x_old[3];
   
	   list<Particle *>::iterator it = particles.begin();
	   double sum=0;
   
	   while (it!=particles.end())
	   {
		   Particle *part = *it;
		   double dt_rem = dt;

		   /*while there is delta t left in particle push*/
		   while(dt_rem>0)
		   {
			   //save old position
			   vec::copy(x_old,part->pos);

			   //new position, x = x + v*dt
			   vec::mult_and_add(part->pos,part->pos,part->vel,dt_rem);
		   
			   //did particle hit the wall?
			   double rr = part->pos[0]*part->pos[0]+part->pos[1]*part->pos[1];
			   if (rr>=R*R)
			   {
				   /*compute intersection point*/
				   double r = sqrt(rr);
				   double r0 = sqrt(x_old[0]*x_old[0]+x_old[1]*x_old[1]);
				   double t = (R-r0)/(r-r0);
       
				   //push particle to the wall
				   vec::mult_and_add(part->pos,x_old,part->vel,t*dt_rem);

				   //update remaining delta_t
				   dt_rem -= t*dt_rem;
           
				   //surface normal vector
				   double normal[3];
				   normal[0] = 0-part->pos[0];
				   normal[1] = 0-part->pos[1];
				   normal[2] = 0;
           
				   vec::unit(normal,normal);
          
				   // get tangents
				   double tang1[3], tang2[3];
				   tangentsFromNormal(normal,tang1,tang2);

				   //cosine emission
					double dir[3];
					lambertianVector(dir, normal, tang1, tang2);
					for (int i=0;i<3;i++)
						part->vel[i]=dir[i]*vdrift;
			   } /*dt_rem*/
			   else //no wall impct
			   {
				    //particle left through inlet or oulet?
				   if (part->pos[2]<=0 || part->pos[2]>=L)
				   {
					   if (part->pos[2]<0)
						   inlet++;
					   else
						   outlet++;
					   //kill particle
					   delete part;
					   it = particles.erase(it);
					   break;		/*break out of dt_rem loop*/
				   }
				   else    //particle inside the cylinder
				   {
					   dt_rem = 0;
					   it++;	/*go to next particle*/
				   }
			   }
		   } ///while dt_rem		   
       }	//particle loop
	}

	/*generates two tangent vectors, assumes normalized input*/
	void tangentsFromNormal(double n[3], double tang1[3], double tang2[3])
	{
		//get maximum direction of the normal
		int max_i=0;
		double max = abs(n[0]);
		for (int i=1;i<3;i++) if (abs(n[i])>max) {max=abs(n[i]);max_i=i;}

		//create a test vector that is not parallel with norm
		double test[3] = {0,0,0};
		if (max_i<2) test[max_i+1]=1.0; else test[0]=1.0;

		//cross the two vectors, this will give the first tangent -perpendicular to both
		vec::cross(tang1,n,test);	/*assuming n is already normalized*/
	
		//and the second vector//
		vec::cross(tang2,n,tang1);
	}

	/*returns direction from cosine distribution about normal given two tangents*/
	void lambertianVector(double dir[3], double normal[3], double tang1[3], double tang2[3])
	{
		//sample angle from cosine distribution
		double sin_theta = sqrt(rnd());	 //sin_theta = [0,1)
		double cos_theta = sqrt(1-sin_theta*sin_theta);
		
		//random rotation about surface normal
		double phi = 2*PI*rnd();
	
		double vn[3],vt1[3],vt2[3];
		vec::mult(vn,normal,cos_theta);
		vec::mult(vt1,tang1,sin_theta*cos(phi));
		vec::mult(vt2,tang2,sin_theta*sin(phi));

		//add components
		for (int i=0;i<3;i++)
			dir[i] = vn[i]+vt1[i]+vt2[i];
	}
};			



/*********   MAIN *********/
int main(int nargs, char*args[])
{
  int num_threads = thread::hardware_concurrency();

  cout<<"Number of supported concurrent threads: "<<num_threads<<endl;

  double L=1;
  double part_totf = 1e7;	//double so we can use scientific notation

  //get parameters from command line if specified
  if (nargs>1) L=atof(args[1]);
  if (nargs>2) num_threads=atoi(args[2]);
  if (nargs>3) part_totf=atof(args[3]);

  //show active parameters
  cout<<"L = "<<L<<"\tnum_threads = "<<num_threads<<"\tpart_tot = "<<part_totf/1e6<<"mil"<<endl;
  
  long part_tot = (long)part_totf;  
  int num_ts = 1000000;	//number of time steps
  long part_per_thread = (long)(part_tot/((double) num_threads)+0.5);

  long np_per_ts = (long)(part_per_thread/(double)num_ts+0.5);
  if (np_per_ts<1) np_per_ts=1;	//make sure we get at least 1
    
  clock_t start = clock();	//start timer

  list<ConductanceXYZ> sims;	//storage for our simulations

  //create num_threads simulations
  for (int i=0;i<num_threads;i++)
	  sims.emplace_back(L,np_per_ts, num_ts);

  //wait for threads to finish
  bool finished;
  do
  {
	   long tot_gen = 0;
	   long tot_out = 0;
	   long tot_np = 0;
	   int ts = num_ts+1;
	   finished = true;
	   //combine counts across threads and also check for completion
	   for (ConductanceXYZ &sim:sims) {
			tot_gen += sim.inlet + sim.outlet;
			tot_out += sim.outlet;
			tot_np += sim.particles.size();
			if (sim.ts<ts) ts=sim.ts;
			finished &= sim.finished;	//boolean and, any false will clear it
	   }
	   double K=0;
	   if (tot_gen>0) K=tot_out/(double)(tot_gen); 
  	   cout<<"it="<<ts<<"\t np="<<setprecision(4)<<tot_np<<"\t K="<<K<<endl;

	   //sleep for 50 milliseconds
	   this_thread::sleep_for(chrono::milliseconds(50));
	} while (!finished);

  clock_t end = clock();
  cout<<"Simulation took "<<(end-start)/((double)CLOCKS_PER_SEC)<<" seconds"<<endl;
  return 0;
}