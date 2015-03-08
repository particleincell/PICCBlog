/*------ Main.java -----------*/
import java.util.*;
 
public class Main 
{
	public static void main(String[] args) 
	{
		ThreadGroup tg = new ThreadGroup("main");
		int np = Runtime.getRuntime().availableProcessors();
		int i, ns=24;
 
		List<Simulation> sims = new ArrayList<Simulation>();
 
		long start = System.nanoTime();
 
		for (i=0;i<ns;i++)
			sims.add(new Simulation(50000,"PI"+i, tg));
 
		i=0;
		while (i<sims.size())
		{
			/*do we have available CPUs?*/
			if (tg.activeCount()<np)
			{
				Simulation sim = sims.get(i);
				sim.start();
				i++;
			} else
				try {Thread.sleep(100);} /*wait 0.1 second before checking again*/
					catch (InterruptedException e) {e.printStackTrace();}
		}
 
 
		/*wait for threads to finish*/
		while(tg.activeCount()>0)
		{
			try {Thread.sleep(100);} 
				catch (InterruptedException e) {e.printStackTrace();}
		}
 
		/*sum up errors*/
		double sum=0;
		for (i=0;i<sims.size();i++)
		{
			Simulation sim = sims.get(i);
			sum+=sim.getError();
		}
 
		long end  = System.nanoTime();
 
		System.out.printf("Average error is %g\n", sum/sims.size());
		System.out.printf("Simulation took %.2g seconds\n", (double)(end-start)/1e9);
	}
}