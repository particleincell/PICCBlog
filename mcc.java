/* a simple demonstration program of the Monte Carlo Collision (MCC) method
 * for calculation of charge exchange (CEX) collisions as described in
 * http://www.particleincell.com/2011/mcc/
 * 
 * developed by Particle In Cell Consulting, particleincell.com, 
 * info@particleincell.com
 */
package mcc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.LinkedList;

/*Main class, creates and launches a new simulation*/
public class mcc
{
	public static void main(String[] args) 
	{
		new MCC_Example().Run();
	}
}

/*the actual simulation is performed by this class*/
class MCC_Example
{
	class Particle 
	{
		double x,y;
		double u,v;
		boolean cex;
		final double q_m=1.602e-19/m_XE;	/* q/m, m=131.3amu (Xenon)*/
		
		/*creates a new particle*/
		Particle ()
		{
			x = 0;
			y = rT*(-1+2*Math.random());		/*scale by thruster radius*/
			u = 29000;
			v = -2000 + Math.random()*4000;
			cex = false;
		}
		
		/*updates velocity*/
		void accelerate(double dt)
		{
			/*fixed electric field*/
			double Ex = 0, Ey =0;
			
			/*only apply within 4 thruster radii of exit*/
			if (x<4*rT)
			{
				double E0=50;
				double theta = 5*(Math.PI/180);	/*in radians*/
				Ex = -Math.sin(theta)*E0;
				Ey = Math.cos(theta)*E0;
				if (y<0) Ey=-Ey;		/*flip below the z axis*/
			}
			
			/*standard Lorentz force push without a mag. field*/
			u+=q_m*Ex*dt;		
			v+=q_m*Ey*dt;
		}
		
		/*updates position*/
		void move(double dt)
		{
			x+=u*dt;
			y+=v*dt;
		}
		
		/*uses MCC to see if the particle collided, and if so, calls CEX*/
		void collideMCC(double dt)
		{
			/*get neutral density at particle location*/
			double nn = NeutralDensity(x, y);
			
			/*get ion velocity*/
			double g = Math.sqrt(u*u + v*v);
			
			/*calculate CEX cross-section for Xenon using Rapp and Francis*/
			double a = -0.8821*Math.log(g)+15.1262;
			double sigma = a*a*1e-20;
					
			/*calculate collision probability*/
			double P = 1 - Math.exp(-nn*sigma*g*dt);
			
			/*compare to a random number*/
			if (Math.random()<P)
			{
				/*collision occured*/
				performCEX();
			}
		}
		
		/*simple collision handler for CEX*/
		void performCEX()
		{
			u = neut_vth*fmaxw();	/*sample new random thermal velocity*/
			v = neut_vth*fmaxw();
			cex=true;				/*flag the particle as CEX for visualization*/
		}
	};
	
	LinkedList<Particle> particles = new LinkedList<Particle>();
	
	final double dt = 1e-7;		/*time step*/
	final double rT = 0.15;		/*thruster radius*/
	final double m_XE = 2.18e-25;	/*131.3 amu, Xenon mass*/
	final double neut_vth = Math.sqrt(2*1.38e-23/m_XE);
	
	/*samples from the Maxwellian distribution using the method of Birdsall*/
	double fmaxw() {return 2*(Math.random()+Math.random()+Math.random()-1.5);}
	
	void Run() 
	{
		int it;
		
		/*main loop*/
		for (it=1;it<=5000;it++)
		{
			InjectParticles();
			CollideParticles();
			MoveParticles();
			
			if (it%50==0)
				OutputParticles(it);
			
			System.out.printf("it: %d \tnp: %d\n",it, particles.size());
		}
	}
	
	/*adds 25 particles per call*/
	void InjectParticles()
	{
		for (int p=0;p<25;p++)
			particles.add(new Particle());
	}
	
	/*loops over particles and moves them*/
	void MoveParticles()
	{
		Iterator iterator = particles.iterator();
		while(iterator.hasNext())
		{
			Particle part = (Particle)iterator.next();
			
			part.accelerate(dt);
			part.move(dt);
			
			/*make sure the particle is still in domain*/
			if (!InDomain(part))
				iterator.remove();
		}
	}
	
	/*loops through particles and calls the MCC algorith for each*/
	void CollideParticles()
	{
		Iterator iterator = particles.iterator();
		while(iterator.hasNext())
		{
			Particle part = (Particle)iterator.next();
			part.collideMCC(dt);
		}
	}
		
	/*returns false if particle is outside our domain box*/
	boolean InDomain(Particle part)
	{
		if (part.x<-0.5 || part.x>2 || part.y<-1 || part.y>1)
			return false;
		return true;
	}
	
	/*returns neutral density at some z,r point*/
	double NeutralDensity(double z, double r)
	{
		final double n0=1e18;
		final double a = 1/(1-1/(Math.sqrt(2)));
		
		double R = Math.sqrt(r*r + (z+rT)*(z+rT));
		double theta = Math.atan(r/(z+rT));
		
		return n0*a*(1-1/Math.sqrt(1+(rT/R)*(rT/R)))*Math.cos(theta);
	}
	
	/*output file handle*/
	PrintWriter part_file = null;
		
	/*saves particles to a tecplot file*/
	void OutputParticles(int it)
	{
	
		/*first time? if so, open the output file*/
		if (part_file == null)
		{
			try
			{
				part_file = new PrintWriter(new FileWriter("particles.dat"));
			}
			catch (Exception e)
			{
				System.err.println(e.getMessage());
			}
			
			part_file.println("VARIABLES = x y u v cex");
		}
		
		/*output new zone*/
		part_file.println("ZONE I="+particles.size()+" T=IT_"+it);
		Iterator<Particle> iterator = particles.iterator();
		while (iterator.hasNext())
		{
			Particle part = iterator.next();
			part_file.printf("%g %g %g %g %d\n",part.x,part.y,part.u,part.v,
						part.cex?1:0);
			
		}
		part_file.flush();
	}
}
