/*
PIC-C Blog: Surface Interaction Test
https://www.particleincell.com/2017/line-triangle-intersection/
Written by: Lubos Brieda, 12/20/2017
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <list>
#include <map>
#include <chrono>
#include "surf-int.h"
#include "surface.h"

using namespace std;

/*particle*/
struct Particle
{
	double pos[3];
	double vel[3];
	int id;
	bool has_trace;

	Particle(SurfaceElement &element, int id) : id(id)
	{
		//random point in the element
		bool primary = element.randomPos(pos);

		//random velocity direction
		element.lambertianVector(vel,primary);	

		//set some magnitude, normally we would smaple from Maxwellian
		vec::mult(vel,vel,10);

		has_trace = false;
	}
};

class TraceSample
{
public:
	double pos[3];
	double vel[3];
	int ts;
	int el_hit = -1;
	TraceSample (Particle &part, int ts, int el_hit = -1) {this->ts=ts;vec::copy(pos,part.pos);vec::copy(vel,part.vel);}
};

class TraceData 
{
public:
	int id;	/*particle id*/
	bool active;
	vector<TraceSample> samples;

	TraceData(Particle &part) {id = part.id;active=true;part.has_trace=true;}
	
};


class ParticleTrace 
{
public:
	ParticleTrace(){};
	
	void add(Particle &part) {traces.emplace_back(part);};
	void updateSingleParticle(Particle &part, int ts, int el_hit = -1);
	void finish();

	size_t size() {return traces.size();}
protected:
	list<TraceData> traces;	
};


const double PI = acos(-1.0);
std::mt19937 mt_gen(0);		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd() {return rnd_dist(mt_gen);}


/*PROTOTYPES*/
bool LoadUnv(const string &file_name, Surface &surface);
void SurfaceSaveVTK(Surface &surface);
void OutputMesh(int ts, Surface &surface, double *count);
void OutputParticles(vector<Particle> &particles);
void SampleParticle (Particle &part, SurfaceElement &element);


/**************** MAIN **************************/
int main()
{
	/*instantiate volume*/
	Surface surface;
	if (!LoadUnv("box.unv",surface) ||
		!LoadUnv("top.unv",surface)) 
	{
		cerr<<"Failed to open input files"<<endl;
		return (-1);
	}

	surface.init();
	
	double dt0 = 1e-4;
	
	/*particle list*/
	vector<Particle> particles; 

	SurfaceZone *sz = surface.findZone("source");
	if (!sz) {cerr<<"Failed to find zone `source'"<<endl;exit(-1);}

	//memory for traces
	ParticleTrace trace;

	//generate particles
	while (particles.size()<100)
	{
		//pick random element
		int i = (int)(rnd()*sz->elements.size());
		SurfaceElement &ele = sz->elements[i];
		if ((ele.area/sz->area) <= rnd())
		{
			int id = (int) particles.size();
			particles.emplace_back(Particle(ele,id));

			//vector makes a local copy so need to pass data from there
			if (particles.size()<5 || id==58)
				trace.add(particles.back());
		}
	}
	
	int np0 = (int) particles.size();
	int total_hits = 0;

	/*timing*/
	auto start = chrono::high_resolution_clock::now();
	
	/*main loop*/
	int ts;
	for (ts=1;ts<=50000;ts++)
	{
		if (ts%100==0) 
			cout<<"TS: "<<ts<<", np: "<<particles.size()<<endl;

		/*update particle positions*/
		vector<Particle>::iterator it = particles.begin();
		while (it!=particles.end())
		{			
			SurfaceElement *el_min = nullptr;
			Particle &part = *it;
			
			int bounces = 0;
			const int max_bounces = 100;
			double part_dt = dt0;
			bool alive  = true;
		
			/*particle loop*/
			while(part_dt>0  && bounces++<max_bounces)
			{
				double x_old[3];
				vec::copy(x_old,part.pos);

				vec::mult_and_add(part.pos,part.pos,part.vel,part_dt);

				//placeholder for retrieving elements from octree or a similar structure
				vector<SurfaceElement*> elements = surface.getElements(x_old,part.pos);
				
				/*intersect with surface triangles*/
				double t_min = 2;
				bool secondary_min = false;
				el_min=nullptr;

				for (SurfaceElement *el:elements)
				{			
					//obtain parametric position of the intersection point
					double t = el->lineIntersect(x_old,part.pos,t_min);

					//check for intersection with secondary triangle on warped quads
					bool secondary=false;
					if (t>=2.0) {t-=2.0;secondary=true;}

					//small positive offset to avoid self intersection with starting surfaces
					if (t<1e-8) continue;

					/*accept if closer than previous hit or if same but looking at us.*/						
					if (t<t_min || (t==t_min && vec::dot(part.vel,el->normal)<vec::dot(part.vel,el_min->normal))) 
					{
						t_min = t; el_min = el; secondary_min = secondary;										
					}				
				}

				/*if no particle impact*/
				if (el_min==nullptr) 
				{
					part_dt=0;
				}
				/*we have impact!*/
				else				
				{			
					//update surface hit counter
					el_min->num_hits++;
					total_hits++;

					//push to the surface, doing this also for sticking particle to capture trace
					vec::mult_and_add(part.pos,x_old,part.vel,part_dt*t_min);
					part_dt -= part_dt*t_min;	

					//compute post impact velocity
					double vel_mag = vec::mag(part.vel);

					/*flip if we hit on back side, make sure to use correct normal on quads*/							
					double *normal = secondary_min?el_min->sec_normal:el_min->normal;
					if (vec::dot(part.vel,normal)>0)
						vel_mag *=-1.0;						
						
					/*diffuse reflection*/
					double ref_vec[3];
					el_min->lambertianVector(ref_vec, !secondary_min);
					vec::mult(part.vel,ref_vec,vel_mag);
					
					/*make sure the point is not on edge*/
					el_min->pushOffEdge(part.pos,!secondary_min);
						
					//add bounce location to trace data (placed here after possible shift from edges)
					if (part.has_trace)
					{								
						trace.updateSingleParticle(part, ts, el_min->id);
					}
						
					}	//if surface hit
				
			}	/*dt loop*/

			//add bounce location to trace data (placed here after possible shift from edges)
			if (part.has_trace)
			{								
				trace.updateSingleParticle(part, ts);
			}

			//check if still inside the bounding box
			if (!alive || !surface.inBoundingBox(part.pos,0.01)) 
			{
				cout<<"Removing particle "<<part.id<<" at ts = "<<ts<<endl;
				it = particles.erase(it);
			}
			else it++;
		}
	} 
	
	int np = (int)particles.size();
	auto end = chrono::high_resolution_clock::now();
	long long interval = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	cout<<"Time per particle push: "<<(double)interval/(ts*np)<<"ms"<<endl;

	//info about particle leaks
	cout<<"Number of leaked particles: "<<(np0-np)<<endl;
	if ((np0-np)>0)
		cout<<"Leak probability: 1 in "<<setprecision(4)<<(double)total_hits/(np0-np)<<endl;
		
	//save trace data
	trace.finish();

	//save surface
	SurfaceSaveVTK(surface);

	return 0;
}

/*** FUNCTIONS ***/
namespace StringUtils
{
	//removes leading and trailing spaces
	std::string trim(const std::string &s)
	{
		const string delims = " \f\n\r\t\v";
		string str= s.substr( 0, s.find_last_not_of(delims) + 1 );
		if (str.empty()) return str;
		str = str.substr( s.find_first_not_of(delims) );
		return str;
	}

	//splits the string by spaces
	std::vector<std::string> split(const std::string &input)
	{
		std::stringstream ss(trim(input));
		std::string s;
		std::vector<std::string> out;

		while(ss) 
		{
			ss>>s;
			out.emplace_back(s);
		}
		return out;
	}

};


/*loads surface mesh from a .unv file*/
bool LoadUnv(const string &file_name, Surface &surf)
{
	ifstream in(file_name);
	if (!in) {return false;}

	cout<<"Reading "<<file_name<<endl;
	/*surface zone*/
	SurfaceZone *surf_zone = nullptr;

	/*lut for index numbers*/
	map<int,int>node_lut;

	int n_nodes = 0;
	int n0 = (int)surf.nodes.size();
	
	string line;
	/*skip 19 lines*/
	for (int i=0;i<19;i++)
	{
		getline(in,line);
	}

	/*now read in nodes*/
	while (in)
	{
		
		/*read nodes*/
		string line1,line2;
		getline(in, line1);
		getline(in, line2);
		
		stringstream s1(line1);
		int node_id;
		s1>>node_id;

		if (node_id==-1) break;

		stringstream s2(line2);
		double pos[3];
		s2>>pos[0]>>pos[1]>>pos[2];

		/*add offset*/
		surf.addNode(pos[0],pos[1],pos[2]);
		n_nodes++;

		if (node_id!=n_nodes) 
			cerr<<"Non-sequential node numbers"<<endl;
	}
	
	/*skip one line*/
	getline(in, line);

	/*temporary data structure for elements*/
	map<int,SurfaceElement> elements;
	map<int,bool> added;
	vector<string> pieces;
	while (in)
	{
		
		/*read element info line*/
		getline(in, line);

		pieces = StringUtils::split(line);
		if (stoi(pieces[0])==-1) break;

		int type = stoi(pieces[5]);

		SurfaceElement ele;
		ele.id = stoi(pieces[0]);

		if (type == 2)	/*edge segment, ignore this one*/
		{
			/*this one has 3 lines for some reason*/
			getline(in, line);
			getline(in, line);
			ele.num_nodes=2;
		}
		else if (type == 3 || type==4)
		{
			getline(in, line);
			pieces = StringUtils::split(line);
			ele.num_nodes = type;
			for (int i=0;i<type;i++)
				ele.connect[i] = n0+stoi(pieces[i])-1;			
		}			
		else {cerr<<"Unknown element type "<<type<<endl;}

		elements[ele.id]=ele;
		added[ele.id]=false;
		
	}

	/*skip two lines*/
	for (int i=0;i<2;i++)	getline(in,line);

	/*now read in groups*/
	while (in)
	{
		getline(in,line);
		pieces =  StringUtils::split(line);
		if (stoi(pieces[0])==-1) break;
		
		/*8th piece contains number of elements in this group*/
		int n_ele = stoi(pieces[7]);
		string name;
		getline(in, name);
		cout<<"Group: "<<name<<endl;
		SurfaceZone &zone = surf.addZone(StringUtils::trim(name));
		for (int i=0;i<n_ele;i++)
		{
			/*there are 4 entries per group, second one contains the element id*/
			int v1,v2,v3,v4;
			in>>v1>>v2>>v3>>v4;
			SurfaceElement &ele = elements[v2];
			added[v2] = true;
			if (ele.num_nodes==3 || ele.num_nodes==4)
				zone.addElement(ele);
		}
		getline(in,line);	/*read in rest of this line (new line character*/

	}

	/*are there some remaining unasigned elements?*/
	std::map<int, bool>::iterator iter;
	SurfaceZone zone("unassigned");
    for (iter = added.begin();iter!=added.end();iter++)
	{
		/*element that hasn't been added yet*/
		if (iter->second==false)			
		{
			SurfaceElement ele = elements[iter->first];
			if (ele.num_nodes==3 || ele.num_nodes==4)
				zone.addElement(ele);
		}
	}
	if (!zone.elements.empty()) surf.addZone(zone);
     
	cout << "Added " << surf.nodes.size() << " nodes and "<<elements.size()<<" elements"<<endl;	
	in.close();
	return true;
}

/*saves surface mesh*/
void SurfaceSaveVTK(Surface &surf)
{
	ofstream out("surface.vtp");	
	if (!out.is_open()) {cerr<<"Failed to open output file"<<endl;return;}

	out<<setprecision(4);

	out<<"<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">\n";
	out<<"<PolyData>\n";

	/*count total number of elements*/
	int n_eles =0;
	for (SurfaceZone &sz: surf.surface_zones) n_eles += (int)sz.elements.size();
	
	out<<setprecision(3);
	int num_elements = 0;
	for (SurfaceZone &sz:surf.surface_zones) num_elements += sz.elements.size();

	out<<"<Piece NumberOfPoints=\""<<surf.nodes.size()<<"\" NumberOfPolys=\""<<num_elements<<"\">\n";

	/*save points*/
	out<<"<Points>\n";
	out<<"<DataArray Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (unsigned int i=0;i<surf.nodes.size();i++)
	{
		Node &node = surf.nodes[i];
		out<<scientific<<setprecision(5)<<node.pos[0]<<" "<<node.pos[1]<<" "<<node.pos[2]<<"\n";
	}
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	out<<"<Polys>\n";
	out<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (SurfaceZone &sz:surf.surface_zones)
		for (SurfaceElement &ele:sz.elements)
			for (int i=0;i<ele.num_nodes;i++) 
				out<<ele.connect[i]<<" ";
	out<<"\n</DataArray>\n";
	
	int offset = 0;		
	out<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (SurfaceZone &sz:surf.surface_zones)
		for (SurfaceElement &ele:sz.elements)
		{
			offset+=ele.num_nodes;	
			out<<offset<<" ";
		}
	out<<"\n</DataArray>\n";
	out<<"</Polys>\n";
	
		
	/*cell data*/
	out<<"<CellData Normals=\"normals\" Scalars=\"num_hits\">\n";
		
	//cell ids
	out<<"<DataArray Name=\"element_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	for (SurfaceZone &sz:surf.surface_zones)
		for (SurfaceElement &el:sz.elements)
			out<<el.id<<" ";
	out<<"\n</DataArray>\n";
		
	/*output normals*/
	out<<"<DataArray Name=\"normals\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (SurfaceZone &sz:surf.surface_zones)
		for (SurfaceElement &el:sz.elements)
			out<<el.normal[0]<<" " <<el.normal[1]<<" "<<el.normal[2]<<"\n";
	out<<"</DataArray>\n";

	out<<"<DataArray Name=\"zone_num\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	int zone_num = 0;
	for (SurfaceZone &sz:surf.surface_zones)
	{
		for (SurfaceElement &el:sz.elements)
			out<<zone_num<<" ";
		zone_num++;
	}
	out<<"\n</DataArray>\n";

	/*surface hit data*/

	out<<"<DataArray Name=\"hit rate\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";		
	for (SurfaceZone &sz:surf.surface_zones)
		for (SurfaceElement &el:sz.elements)
			out<<el.num_hits/el.area<<" ";			
		
	out<<"\n";
	out<<"</DataArray>\n";
 out<<"</CellData>\n";
 out<<"</Piece>\n";
 out<<"</PolyData>\n";
 out<<"</VTKFile>\n";
 out.close();
}



//particle trace

void ParticleTrace::updateSingleParticle(Particle &part, int ts, int el_hit)
{
	/*update all traces*/
	for (TraceData &td:traces)
		if (td.id==part.id)
		{
			if (!td.active) return;	
			td.samples.emplace_back(part, ts, el_hit);			
			return;
		}
}

double format(double in)
{
	if (fabs(in)<1e-15) return 0.0;
	return in;

}

/*outputs trace data to a file*/
void ParticleTrace::finish()
{
	ofstream out("trace.vtp");
	if (!out.is_open()) return;

	long np = 0;
	for (TraceData &td:traces)		
		np += (int)td.samples.size();

	out<<"<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">\n";
	out<<"<PolyData>\n";
	out<<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfLines=\""<<traces.size()<<"\">\n";

	out<<setprecision(4)<<scientific;
	
	out<<"<Points>\n";
	out<<"<DataArray Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (TraceData &td:traces)
		for (TraceSample &s:td.samples)
			out<<format(s.pos[0])<<" "<<format(s.pos[1])<<" "<<format(s.pos[2])<<"\n";
	out<<"</DataArray>\n";
	out<<"</Points>\n";
	
	/*point data*/
	out<<"<PointData Scalars=\"time\" Vectors=\"vel\">\n";
	/*velocity*/
	out<<"<DataArray Name=\"vel\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (TraceData &td:traces)
		for (TraceSample &s:td.samples)
			out<<format(s.vel[0])<<" "<<format(s.vel[1])<<" "<<format(s.vel[2])<<"\n";
	out<<"</DataArray>\n";

	out<<"<DataArray Name=\"time step\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	for (TraceData &td:traces)
		for (TraceSample &s:td.samples)
			out<<s.ts<<" ";
	out<<"</DataArray>\n";

	out<<"<DataArray Name=\"trace_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	for (TraceData &td:traces)
		for (TraceSample &s:td.samples)
			out<<td.id<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"element_hit\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	for (TraceData &td:traces)
		for (TraceSample &s:td.samples)
			out<<s.el_hit<<" ";
	out<<"\n</DataArray>\n";
	out<<"</PointData>\n";
		
	out<<"<Lines>\n";
	out<<"<DataArray type=\"Int32\" Name=\"connectivity\">\n";
	long n=0;
	for (TraceData &td:traces)
		for (int i=0;i<(int)td.samples.size();i++) out<<n++<<" ";
	out<<"\n</DataArray>\n";
	out<<"<DataArray type=\"Int32\" Name=\"offsets\">\n";
	n=0;
	for (TraceData &td:traces)
	{	n+=(int)td.samples.size();
		out<<n<<" ";	
	}
	out<<"\n</DataArray>\n";
	out<<"</Lines>\n";

	out<<"</Piece>\n";

 out<<"</PolyData>\n";
 out<<"</VTKFile>\n";

	out.close();
}
