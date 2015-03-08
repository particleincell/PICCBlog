/*
 * ParticleSpeedTest
 * 
 * Code to test several different particle data structures
 * See: http://www.particleincell.com/2012/particle-data-structure/
 */

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

public class ParticleSpeedTest 
{
	static class Field
	{
		double data[][][];
		
		int ni,nj,nk;		/*number of nodes*/
		double dh;			/*node spacing*/
		
		double x0[],xd[];	/*origin and diagonal corner*/
		
		/*constructor*/
		Field(int ni, int nj, int nk, double dh)
		{
			/*allocate memory*/
			data = new double[ni][][];
			for (int i=0;i<ni;i++)
			{
				data[i] = new double[nj][];
				for (int j=0;j<nj;j++)
					data[i][j] = new double[nk];
			}
			
			/*save dimensions*/
			this.ni = ni;
			this.nj = nj;
			this.nk = nk;
			this.dh = dh;
			
			/*set mesh corners, defaulting to (0,0,0) origin*/
			x0 = new double[3];
			xd = new double[3];
			x0[0]=0;x0[1]=0;x0[2]=0;
			
			/*diagonal corner, (ni-1) is the number of cells*/
			xd[0]=x0[0]+(ni-1)*dh;
			xd[1]=x0[1]+(nj-1)*dh;
			xd[2]=x0[2]+(nk-1)*dh;
		}
		
		/*returns true if point is in domain, x=[x0,xd) (right edge exclusive)*/
		boolean inDomain(double x[])
		{
			for (int dim=0;dim<3;dim++)
				if (x[dim]<x0[dim] || x[dim]>=xd[dim])
					return false;
			
			return true;
		}
		
		/*converts physical coordinate to a logical one*/
		double[] XtoL (double x[])
		{
			double lc[] = new double[3];
			
			for (int dim=0;dim<3;dim++)
				lc[dim]=(x[dim]-x0[dim])/dh;
			
			return lc;
		}
		
		/*clears all values*/
		void clear()
		{
			for (int i=0;i<ni;i++)
				for (int j=0;j<nj;j++)
					Arrays.fill(data[i][j], 0);
		}
		
		/*scatters value onto logical coordinate lc[]*/
		void scatter(double val, double lc[])
		{
			int i = (int)lc[0];
			int j = (int)lc[1];
			int k = (int)lc[2];
			double di = lc[0]-i;
			double dj = lc[1]-j;
			double dk = lc[2]-k;
			
			/*standard linear weighing*/
			data[i][j][k] +=  (1.0-di)*(1.0-dj)*(1.0-dk);
			data[i][j][k+1] +=  (1.0-di)*(1.0-dj)*dk;
			data[i][j+1][k] += (1.0-di)*(dj)*(1.0-dk);
			data[i][j+1][k+1] += (1.0-di)*(dj)*dk;
		
			data[i+1][j][k] += di*(1.0-dj)*(1.0-dk);
			data[i+1][j][k+1] += di*(1.0-dj)*dk;
			data[i+1][j+1][k] += di*dj*(1.0-dk);
			data[i+1][j+1][k+1] += di*dj*dk;
		}
		
		/*outputs the field in the Tecplot ASCII format*/
		void save(String file_name)
		{
			PrintWriter pw = null;
			
			try
			{
				pw = new PrintWriter(new FileWriter(file_name));
			}
			catch (IOException e)
			{
				System.err.println("Failed to open output file "+file_name);
			}
			
			/*write header*/
			pw.println("VARIABLES = x y z value");
			pw.printf("ZONE I=%d J=%d K=%d\n",ni,nj,nk);
			
			/*write data assuming (0,0,0) origin*/
			for (int k=0;k<nk;k++)
				for (int j=0;j<nj;j++)
					for (int i=0;i<ni;i++)
						pw.printf("%g %g %g %g\n", x0[0]+i*dh, x0[1]+j*dh, 
												   x0[2]+k*dh, data[i][j][k]);
				
			/*close file*/
			pw.close();
		}
	}
	
	static class Particle
	{
		double x[] = new double[3];
		double v[] = new double[3];
	}
	
	static class ParticleBlockList extends LinkedList<Particle>
	{
		protected LinkedList<ParticleBlock> part_block = new LinkedList<ParticleBlock>();
		protected static final int block_size = 5000;			/* maximum number of particles in a block*/
		protected int np=0;
		
		@Override
		public int size() {return np;}
		
		/** Adds particle specified by x and v*/
		
		@Override
		public boolean add(Particle part)
		{
			/*add block if one doesn't exist already*/
			if (part_block.size()==0)
				part_block.add(new ParticleBlock());

			ListIterator<ParticleBlock> iterator = (ListIterator<ParticleBlock>)part_block.iterator();
			ParticleBlock block = iterator.next();
			
			/*loop through blocks looking for an empty slot*/
			while (true)
			{
				if (block.count<block_size)
				{
					block.particle[block.count] = part;
					block.count++;
					break;
				}
				else		/*no room in this block*/
				{
					/*go to next block*/
					if (iterator.hasNext())
						block = iterator.next();
					else /*add a new block if we have reached the end*/
					{
						block = new ParticleBlock();
						iterator.add(block);
					}
				}
			} /*while*/

			/*increment counter*/
			np++;

			return true;
	}

		/** Particle block*/
		public static class ParticleBlock
		{
			int count = 0;										/* number of particles*/
			Particle particle[] = new Particle[block_size];		/* particle array*/
		};

		/** Iterator for moving through particles*/
		public class ParticleIterator implements Iterator<Particle>
		{
			protected Particle part;		
			protected ParticleBlock curr_block;		/* index of last particle returned by NextParticle*/
			protected int curr_index;				/* block index for curr_block*/
			protected int curr_block_index;			/* block index for curr_block*/
			protected boolean repeat;

			ParticleIterator()
			{
				if (part_block.size()>0)
					curr_block = part_block.getFirst();
				else
					curr_block = null;

				curr_index = -1;
			}
		
			@Override
			public boolean hasNext() 
			{
				if (curr_block==null) return false;

				if (curr_index+1<curr_block.count ||
					curr_block_index<part_block.size()-1) return true;
				else return false;
			}

			@Override
			/** returns the next particle*/
			public Particle next() 
			{
				assert(curr_block!=null);

				if (repeat==false)
				{
					/*advance next pointer*/
					curr_index++;

					if (curr_index>=curr_block.count)
					{
						curr_block_index++;
						curr_index=0;

						if (curr_block_index<part_block.size())
							curr_block = part_block.get(curr_block_index);
						else
							curr_block = null;	/*we reached the end*/
					}
				}

				repeat=false;

				part = curr_block.particle[curr_index];
				return part;
			}

			@Override
			/** removes current particle*/
			public void remove() 
			{
				/*may want to do a clone() here*/
				if (curr_block.count>1)
					curr_block.particle[curr_index] = curr_block.particle[curr_block.count-1];

				curr_block.count--;
				repeat = true;

				/*did we remove the sole particle from this block?*/
				if (curr_block.count==0)
				{
					if (curr_block_index==part_block.size()-1 &&
						curr_block_index>0	)			/*last block*/
						curr_block_index--;

					if (part_block.size()>1)
						part_block.remove(curr_block);

					curr_block = part_block.get(curr_block_index);
				}

				/*did we remove the last particle?*/
				if (curr_index>=curr_block.count)
				{
					repeat=false;
				}

				np--;	
			}
		} /*iterator*/
		
		/** returns new iterator*/
		public ParticleIterator iterator() {return new ParticleIterator();}
	}
	
	/*main loop*/
	public static void main(String args[])
	{
		/*allocate density field*/
		Field den = new Field(30,30,30,0.05);
		
	//	ArrayList<Particle> particles = new ArrayList<Particle>();
		LinkedList<Particle> particles = new LinkedList<Particle>();
	//	ParticleBlockList particles = new ParticleBlockList();
		
		final double dt = 1e-5;
		final double spwt = 1000;	/*specific weight (scaling factor)*/
		
		final int max_it=500;
		int np_count[] = new int[max_it];
		
		long start = System.nanoTime();
		
		/*iterate*/
		for (int it=0;it<max_it;it++)
		{
			/*clear density*/
			den.clear();
			
			/*sample particles*/
			SampleParticles(particles,1000);
			
			
			/*loop over particles*/
			Iterator<Particle> iterator = particles.iterator();
			
			while (iterator.hasNext())
			{
				Particle particle = iterator.next();
		
				/*update position*/
				particle.x[0] += particle.v[0]*dt;
				particle.x[1] += particle.v[1]*dt;
				particle.x[2] += particle.v[2]*dt;
				
				/*deposit to grid if in domain, remove otherwise*/
				if (den.inDomain(particle.x))
				{
					double lc[] = den.XtoL(particle.x);
					den.scatter(spwt, lc);
				}
				else	/*remove particle*/
					iterator.remove();
			}
	
			np_count[it]=particles.size();
			if (it%25==0) System.out.println("IT:"+it+" NP:" + particles.size());
		}
		
		long end=System.nanoTime();
		System.out.printf("Simulation took %.4g seconds\n",(end-start)/1e9);
		
		/*save field*/
		den.save("density.dat");
		
		/*save particle count*/
		save1D("np_count.dat",np_count);
	}
	
	static void save1D(String file_name, int data[])
	{
		PrintWriter pw = null;
			
		try	{
			pw = new PrintWriter(new FileWriter(file_name));
		}
		catch (IOException e)
		{
			System.err.println("Failed to open output file "+file_name);
		}
			
		/*write header*/
		pw.println("VARIABLES = i val");
		pw.printf("ZONE I=%d \n",data.length);
		for (int i=0;i<data.length;i++) 
			pw.printf("%d %d\n", i,data[i]);
		pw.close();		
	}
	
	/*simple source to sample the specified number of particles*/
	static void SampleParticles(List<Particle>particles, int count)
	{
		for (int i=0;i<count;i++)
		{
			Particle new_part = MaxwPart();
			particles.add(new_part);
		}
	}
	
	/*returns a particle with the maxwellian velocity*/
	static Particle MaxwPart()
	{
		Particle part = new Particle();
		
		part.x[0] = .75;
		part.x[1] = .75;
		part.x[2] = .75;
		
		final double v_th = 600;	/*some thermal velocity*/
		
		/*Birdsall's method for sampling Maxwellian with M=3*/
		double v_maxw = v_th*(Math.random()+Math.random()+Math.random() - 1.5);
		v_maxw = v_th;
		
		/*method 1*/
		double n[] = new double[3];
		n[0] = -1.0 + 2*Math.random();
		n[1] = -1.0 + 2*Math.random();
		n[2] = -1.0 + 2*Math.random();
		n[2]=0;
		
		double n_mag = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
		n[0] /= n_mag;
		n[1] /= n_mag;
		n[2] /= n_mag;
		
		/*method 2*/
		double theta = 2*Math.PI*Math.random();
		double phi=Math.acos(-1+2*Math.random());
		
		double R= -1.0+2*Math.random();
		double rs=Math.sqrt(1-R*R);
		
		n[0] = Math.cos(theta)*Math.sin(phi);
		n[1] = Math.sin(theta)*Math.sin(phi);
		n[2] = Math.cos(phi);
	
		n[0] = Math.cos(theta)*rs;
		n[1] = Math.sin(theta)*rs;
		n[2] = R;

		part.v[0] = n[0]*v_maxw;
		part.v[1] = n[1]*v_maxw;
		part.v[2] = n[2]*v_maxw;	/*drift speed*/
		
		return part;
	}

}
