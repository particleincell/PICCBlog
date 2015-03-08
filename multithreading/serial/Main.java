/*------ Main.java -----------*/
public class Main 
{
 
	public static void main(String[] args) 
	{
		double sum = 0;
		int ns=24;   /*number of computations*/
 
                /*sample 50,000 x 1,000 points*/
		Simulation sim = new Simulation(50000);
 
		for (int i=0;i<ns;i++)
		{
			sim.computePI();
			sum += sim.getError();
		}
 
		System.out.printf("Average error is %g\n", sum/ns);
	}
}
