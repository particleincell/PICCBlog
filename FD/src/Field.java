/*very very basic 2D mesh */
public class Field {
	protected double[][] data;
	protected int nx,ny;
	
	public Field(int nx, int ny)
	{
		this.nx = nx;
		this.ny = ny;
		
		data = new double[nx][ny];
	}
	
	void set(double val)
	{
		for (int j=0;j<ny;j++)
			for (int i=0;i<nx;i++) data[i][j]=val;
	}
	
	double[][] getData() {return data;}
	double at(int i, int j) {return data[i][j];}
}
