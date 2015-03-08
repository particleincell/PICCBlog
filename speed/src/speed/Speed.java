package speed;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;

/* ************************************************************
 * Code to test the performance of a finite difference 
 * Poisson solver based on memory access ordering
 * 
 * For more information see:
 * http://www.particleincell.com/2012/memory-code-optimization/
 * 
 * ***********************************************************/
public class Speed 
{
	final static int nn = 50;
	final static int num_it = 500;
	final static int nx=nn,ny=nn,nz=nn;
	
	public static void main(String[] args) 
	{
		int i,j,k;
		int it;
		long start, end;
		
		/*array*/
		int size = nx*ny*nz;
		double array[] = new double[size];
		
		for (i=0;i<size;i++)
			array[i]=i;
		
		double sum=0;
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
			for (i=0;i<size;i++)
				sum+=array[i];
		end = System.nanoTime();
		System.out.println("Array took "+(1e-9)*(end-start)+" seconds");
		
		/*array list*/
		ArrayList<Double> array_list = new ArrayList<Double>(size);
		
		for (i=0;i<size;i++)
			array_list.add(new Double(i));
		
		sum=0;
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
			for (i=0;i<size;i++)
				sum+=array_list.get(i);
		end = System.nanoTime();
		System.out.println("ArrayList took "+(1e-9)*(end-start)+" seconds");
		
		/*linked list*/
		sum = 0;
		LinkedList<Double> linked_list = new LinkedList<Double>();
		for (i=0;i<size;i++)
			linked_list.add(new Double(i));
		
		start = System.nanoTime();
		for (it=0;it<20;it++)
		{
			for (i=0;i<size;i++)
				sum+=linked_list.get(i);
		}
		end = System.nanoTime();
		System.out.println("LinkedList took "+(1e-9)*(num_it/it)*(end-start)+" seconds");
		
		
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
			for (i=0;i<size;i++)
				array[i]=i;
		end = System.nanoTime();
		System.out.println("Array took "+(1e-9)*(end-start)+" seconds");
		
		/*data structures for 3D and 1D approaches*/
		double x3[][][] = allocate3D();
		double x1[] = allocate1D();
			
		/*case 1, 3D k->j->i*/
		initData3D(x3);	
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
		{
			for (k=1;k<nz-1;k++)
				for (j=1;j<ny-1;j++)
					for (i=1;i<nx-1;i++)
					{
						x3[i][j][k]=(1/6.0)*(x3[i-1][j][k]+x3[i+1][j][k]+
											 x3[i][j-1][k]+x3[i][j+1][k]+
											 x3[i][j][k-1]+x3[i][j][k+1]);
					}
		}
		end = System.nanoTime();
		System.out.println("Case 1, 3D k->j->i took "+(1e-9)*(end-start)+" seconds");
		output3D("case1.dat",x3);
		
		/*case 2, 3D i->j->k*/
		initData3D(x3);	
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
		{
			for (i=1;i<nx-1;i++)
				for (j=1;j<ny-1;j++)
					for (k=1;k<nz-1;k++)
					{
						x3[i][j][k]=(1/6.0)*(x3[i-1][j][k]+x3[i+1][j][k]+
											 x3[i][j-1][k]+x3[i][j+1][k]+
											 x3[i][j][k-1]+x3[i][j][k+1]);
					}
		}
		end = System.nanoTime();
		System.out.println("Case 2, 3D i->j->k took "+(1e-9)*(end-start)+" seconds");
		output3D("case2.dat",x3);
	
		/*case 3, 1D flat array*/
		initData1D(x1);	
		/*precompute node offsets*/
		int node_offset[] = nodeOffsets();
		
		start = System.nanoTime();
		for (it=0;it<num_it;it++)
		{
			for (i=1;i<nx-1;i++)
				for (j=1;j<ny-1;j++)
					for (k=1;k<nz-1;k++)
					{
						int u=IJKtoU(i,j,k);
						sum=0;
						for (int t=0;t<6;t++) sum+=x1[u+node_offset[t]];
						
						x1[u]=(1/6.0)*sum;
					}
		}
		end = System.nanoTime();
		System.out.println("Case 3, flat 1D took "+(1e-9)*(end-start)+" seconds");
		output1D("case3.dat",x1);
	}
	
	/**allocates 3D nn*nn*nn array*/
	static double[][][] allocate3D()
	{
		double x[][][]=new double[nx][][];
	
		for (int i=0;i<nx;i++)
		{
			x[i] = new double[ny][];
			for (int j=0;j<ny;j++)
				x[i][j] = new double[nz];	
		}
		return x;		
	}
	
	/**resets data, assumes uniform 3D nn*nn*nn mesh*/
	static void initData3D(double x[][][])
	{
		int i,j,k;
		
		/*set everything to zero*/
		for (i=0;i<nx;i++)
			for (j=0;j<ny;j++)
				for (k=0;k<nz;k++)
					x[i][j][k]=0;
		
		/*set default value of 100 on x=0 plane*/
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
				x[0][j][k]=100;
	}

	/**allocates 1D nn*nn*nn array*/
	static double[] allocate1D()
	{
		/*allocate memory structure*/
		return new double[nn*nn*nn];
	}
	
	/**resets data, assumes uniform 3D nn*nn*nn mesh*/
	static void initData1D(double x[])
	{
		int i,j,k;
		
		/*set everything to zero*/
		for (i=0;i<nx;i++)
			for (j=0;j<ny;j++)
				for (k=0;k<nz;k++)
				{
					x[IJKtoU(i,j,k)]=0;
				}
		
		/*set default value of 100 on x=0 plane*/
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			{
				x[IJKtoU(0,j,k)]=100;
			}
	}
	
	/**flattens i,j,k index, u = i*(ny*nz)+j*(nz)+k*/
	static int IJKtoU(int i, int j, int k)
	{
		return i*(ny*nz) + j*nz + k;
	}
	
	/**returns node offsets for a standard finite difference stencil*/
	static int[] nodeOffsets()
	{
		int node_offsets[] = new int[6];
		node_offsets[0]=IJKtoU(-1,0,0);
		node_offsets[1]=IJKtoU(+1,0,0);
		node_offsets[2]=IJKtoU(0,-1,0);
		node_offsets[3]=IJKtoU(0,+1,0);
		node_offsets[4]=IJKtoU(0,0,-1);
		node_offsets[5]=IJKtoU(0,0,+1);
		return node_offsets;
	}
	
	/**saves 3D mesh in the Tecplot format*/
	static void output3D(String file_name, double x3[][][])
	{
		PrintWriter pw = null;
		try{
			pw = new PrintWriter(new FileWriter(file_name));
		}
		catch (Exception e)
		{
			System.err.println("Failed to open output file "+file_name);
		}
		
		pw.println("VARIABLES = i j k X");
		pw.printf("ZONE I=%d J=%d K=%d\n",nx,ny,nz);
		
		for (int i=0;i<nx;i++)
			for (int j=0;j<ny;j++)
				for (int k=0;k<nz;k++)
					pw.printf("%d %d %d %g\n",i,j,k,x3[i][j][k]);
		pw.close();
	}

	/**saves flat 1D mesh in the Tecplot format*/
	static void output1D(String file_name, double x1[])
	{
		PrintWriter pw = null;
		try{
			pw = new PrintWriter(new FileWriter(file_name));
		}
		catch (Exception e)
		{
			System.err.println("Failed to open output file "+file_name);
		}
		
		pw.println("VARIABLES = i j k X");
		pw.printf("ZONE I=%d J=%d K=%d\n",nx,ny,nz);
		
		for (int i=0;i<nx;i++)
			for (int j=0;j<ny;j++)
				for (int k=0;k<nz;k++)
					pw.printf("%d %d %d %g\n",i,j,k,x1[IJKtoU(i,j,k)]);
		pw.close();		
	}
}
