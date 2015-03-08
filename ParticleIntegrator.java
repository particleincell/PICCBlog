package particleintegrator;

import java.io.*;

/* Simple particle integrator to show how to push particles
 * in particle in cell plasma simulation codes
 * 
 * Derivation and overview is located at:
 * http://www.particleincell.com/2011/vxb-rotation/
 * 
 * Written by Lubos Brieda, lubos.brieda@particleincell.com
 * 
 * Visit http://www.particleincell.com/blog for articles on 
 * scientific computing and plasma modeling
 * 
 * Don't forget to reference particleincell.com if you find this work useful
 */
public class ParticleIntegrator 
{
    public static void main(String[] args) 
    {
	double E[],B[];
	
	double dt = 3e-11;
	
	int it;
	PrintWriter pw=null;
	
	/*open output file*/
	try
	{
	    FileOutputStream out = new FileOutputStream("trace.txt");
	    pw = new PrintWriter(out);
	}
	catch(Exception e)
	{
	    System.err.println(e.toString());
	}
	pw.println("it time x y z u v w");
	
	/*sample particle*/
	Particle part = SampleParticle();
	
	/*push velocity back in time by 1/2 dt*/
	E=EvalE(part.x);
	B=EvalB(part.x);
	UpdateVelocity(part,E,B,-0.5*dt);
	
	long start_time = System.nanoTime();
	
	for (it=0;it<1000;it++)
	{
	    E=EvalE(part.x);
	    B=EvalB(part.x);
	    UpdateVelocity(part,E,B,dt);
	    PushParticle(part,dt);
	    
	    /*output data every 2 time steps*/
	    if (it%2==0) 
	    	pw.printf("%d %g %g %g %g %g %g %g\n",it,it*dt,part.x[0],part.x[1],part.x[2],
						    part.v[0],part.v[1],part.v[2]);	    
	}
	
	long end_time = System.nanoTime();
	
	pw.close();
	System.out.printf("Finished after %d time steps in %g seconds\n",it,
		    (double)(end_time - start_time)*1e-9);
    }
    
    /*updates velocity, acts as a method selector*/
    static void UpdateVelocity(Particle part, double E[], double B[], double dt)
    {
//	UpdateVelocityForward(part, E, B, dt);
//   	UpdateVelocityTajimaImplicit(part, E, B, dt);
//    	UpdateVelocityTajimaExplicit(part, E, B, dt);
    	UpdateVelocityBoris(part, E, B, dt);   	
    }
  
    static class Particle
    {
	double x[] = new double [3];
	double v[] = new double [3];
	final double q = -1.602e-19;	    /*particle charge (electron)*/
	final double m = 9.109e-31;	    /*particle mass (electron)*/
    }

    /*evaluates electric field at particle position, constant in this example*/
    static Particle SampleParticle()
    {
	Particle part = new Particle();
	
	part.x[0] = 0;
	part.x[1] = 0;
	part.x[2] = 0;
	
	part.v[0] = 0;
	part.v[1] = 1e5;
	part.v[2] = 0;
	
	/*sample B to compute larmor radius, this is hardcoded for uniform B=Bz*/
	double B[] = EvalB(part.x);
	double rL = part.m*part.v[1]/(Math.abs(part.q)*B[2]);
	part.x[0] = rL;
	
	System.out.printf("Larmor radius is %g\n",rL);
	return part;
    }

    /*evaluates electric field at particle position, constant in this example*/
    static double[] EvalE(double pos[])
    {
	double E[] = new double[3];
	E[0] = 0;
	E[1] = 0;
	E[2] = 0;
	return E;
    }

    /*evaluates electric field at particle position, constant in this example*/
    static double[] EvalB(double pos[])
    {
	double B[] = new double[3];
	B[0] = 0;
	B[1] = 0;
	B[2] = 0.01;
	return B;
    }
    
    static void PushParticle(Particle part, double dt)
    {
	part.x[0] += part.v[0]*dt;
	part.x[1] += part.v[1]*dt;
	part.x[2] += part.v[2]*dt;
    }
    
    static double[] CrossProduct(double v1[], double v2[])
    {
	double r[]=new double[3];
	r[0] = v1[1]*v2[2]-v1[2]*v2[1];
	r[1] = -v1[0]*v2[2]+v1[2]*v2[0];
	r[2] = v1[0]*v2[1]-v1[1]*v2[0];
	return r;
    }
      
    /*updates velocity using forward differencing: INCORRECT!*/
    static void UpdateVelocityForward(Particle part, double E[], double B[], double dt)
    {
	double vxB[] = CrossProduct(part.v,B);
	
	part.v[0] += part.q/part.m*(E[0] + vxB[0])*dt;
	part.v[1] += part.q/part.m*(E[1] + vxB[1])*dt;
	part.v[2] += part.q/part.m*(E[2] + vxB[2])*dt;
    }

    /*updates velocity using the Tajima method, Computational Plasma Physiscs, p.62-64, INCORRECT*/
    static void UpdateVelocityTajimaExplicit(Particle part, double E[], double B[], double dt)
    {
	double vxB[] = CrossProduct(part.v,B);
	double q_over_m = part.q/part.m;
	double k=part.q/part.m*dt;
	double M[][] = new double[3][3];
	double a[] = new double[3];
	
	/*compute rotation matrix*/
	M[0][0]=1;
	M[0][1]=k*B[2];
	M[0][2]=-k*B[1];
	
	M[1][0] = -k*B[2];
	M[1][1] = 1;
	M[1][2] = k*B[0];
	
	M[2][0] = k*B[1];
	M[2][1] = -k*B[0];
	M[2][2] = 1;
	
	/*compute acceleration vector*/
	for (int dim=0;dim<3;dim++)
	    a[dim] = part.v[dim] + part.q/part.m*0.5*dt*E[dim];
	
	for (int dim=0;dim<3;dim++)
	    part.v[dim] = q_over_m*E[dim]*0.5*dt + M[dim][0]*a[0] + M[dim][1]*a[1] + M[dim][2]*a[2];	
    }			

    /*updates velocity using the Tajima method, Computational Plasma Physiscs, p.62-64*/
    static void UpdateVelocityTajimaImplicit(Particle part, double E[], double B[], double dt)
    {
	double vxB[] = CrossProduct(part.v,B);
	double q_over_m = part.q/part.m;
	double k=part.q/part.m*0.5*dt;
	double M1[][] = new double[3][3];
	double M2[][] = new double[3][3];
	double iM1[][] = new double[3][3];
	double vp[] = new double[3];
	double v1[] = new double[3];
	double v2[] = new double[3];
	double a[] = new double[3];
	
	/*compute [I-Rk]*/
	M1[0][0] = 1;
	M1[0][1] = -k*B[2];
	M1[0][2] = k*B[1];
	
	M1[1][0] = k*B[2];
	M1[1][1] = 1;
	M1[1][2] = -k*B[0];
	
	M1[2][0] = -k*B[1];
	M1[2][1] = k*B[0];
	M1[2][2] = 1;
	
	/*calculate determinant*/
	double det = Determinant(M1);
	
	/*compute inverse [I-Rk]^-1*/
	iM1[0][0] = (1+k*k*B[0]*B[0])/det;
	iM1[0][1] = (k*B[2]+k*k*B[0]*B[1])/det;
	iM1[0][2] = (-k*B[1]+k*k*B[0]*B[2])/det;
	
	iM1[1][0] = (-k*B[2]+k*k*B[0]*B[1])/det;
	iM1[1][1] = (1+k*k*B[1]*B[1])/det;
	iM1[1][2] = (k*B[0]+k*k*B[1]*B[2])/det;
	
	iM1[2][0] = (k*B[1]+k*k*B[0]*B[1])/det;
	iM1[2][1] = (-k*B[0]+k*k*B[1]*B[2])/det;
	iM1[2][2] = (1+k*k*B[2]*B[2])/det;
	
	/*compute [I+Rk]*/
	M2[0][0]=1;
	M2[0][1]=k*B[2];
	M2[0][2]=-k*B[1];
	
	M2[1][0] = -k*B[2];
	M2[1][1] = 1;
	M2[1][2] = k*B[0];
	
	M2[2][0] = k*B[1];
	M2[2][1] = -k*B[0];
	M2[2][2] = 1;
	
	/*compute acceleration vector*/
	for (int dim=0;dim<3;dim++)
	    a[dim]=E[dim]*part.q/part.m*dt;
	
	vp = MatrixVectMult(M2, part.v);	
	v1 = MatrixVectMult(iM1, vp);
	v2 = MatrixVectMult(iM1, a);
	part.v = VectVectAdd(v1, v2);
    }			

    /*updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62*/
    static void UpdateVelocityBoris(Particle part, double E[], double B[], double dt)
    {
	double v_minus[] = new double[3];
	double v_prime[] = new double[3];
	double v_plus[] = new double[3];
	
	double t[] = new double[3];
	double s[] = new double[3];
	double t_mag2;
	
	int dim;
	
	/*t vector*/
	for (dim=0;dim<3;dim++)
	    t[dim] = part.q/part.m*B[dim]*0.5*dt;
	
	/*magnitude of t, squared*/
	t_mag2 = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
	
	/*s vector*/
	for (dim=0;dim<3;dim++)
	    s[dim] = 2*t[dim]/(1+t_mag2);

	/*v minus*/
	for (dim=0;dim<3;dim++)
	    v_minus[dim] = part.v[dim] + part.q/part.m*E[dim]*0.5*dt;
	
	/*v prime*/
	double v_minus_cross_t[] = CrossProduct(v_minus, t);
	for (dim=0;dim<3;dim++)
	    v_prime[dim] = v_minus[dim] + v_minus_cross_t[dim];
	
	/*v prime*/
	double v_prime_cross_s[] = CrossProduct(v_prime, s);
	for (dim=0;dim<3;dim++)
	    v_plus[dim] = v_minus[dim] + v_prime_cross_s[dim];
	
	/*v n+1/2*/
	for (dim=0;dim<3;dim++)
	    part.v[dim] = v_plus[dim] + part.q/part.m*E[dim]*0.5*dt;
    }

    /*determinant of a 3x3 matrix*/
    static double Determinant (double a[][])
    {
	return	a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) -
		a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) + 
		a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    }
    
    /*matrix vector multiplication*/
    static double[] MatrixVectMult(double a[][], double x[])
    {
	double b[] = new double[3];
	
	for (int dim=0;dim<3;dim++)
	    b[dim] = a[dim][0]*x[0] + a[dim][1]*x[1] + a[dim][2]*x[2];
	
	return b;
    }
    
    /*addition of two vectors*/
    static double[] VectVectAdd(double a[], double b[])
    {
	double r[] = new double[3];
	
	for (int dim=0;dim<3;dim++)
	    r[dim]=a[dim]+b[dim];
	
	return r;
    }

}
