import java.io.FileOutputStream;
import java.io.PrintStream;

public class Main {

	/*some globals*/
	static Solver solver;

	static int nx,ny;
	static double dx,dy;
	
	static Field phi,rho;
	
	public static void main(String[] args) 
	{
		/*domain parameters*/
		nx = 51; ny = 21;
		dx = 0.01; dy = 0.02;
				
		/*allocate data*/
		phi = new Field(nx,ny);
		rho = new Field(nx,ny);
		
		/*initialize rho*/
		rho.set(1.602e-19*1e0);
		
		/*initialize field solver*/
		solver = new Solver(phi,rho,dx,dy);
		
		
		solver.update();
			
		saveResults();
			
		System.out.println("Done!");
	
	}
	
	static public void saveResults() 
	{
        int i,j;
        
        PrintStream f_p=null;
    
        try
        {
        	FileOutputStream out = new FileOutputStream("phi.dat");
            f_p = new PrintStream( out );
            f_p.println ("VARIABLES = \"x (m)\" \"y (m)\" phi");
        }
        catch (Exception e)
        {
                System.err.println ("Error writing results");
        }
        
        f_p.printf("ZONE I=%d J=%d \n",nx,ny);
        /*iterate over the field*/
        for (j=0;j<ny;j++)
        	for (i=0;i<nx;i++)
        		f_p.printf("%g %g %g\n", i*dx,j*dy,phi.at(i,j));
        	
        
    }
}

