/*potential solver*/
public class Solver {

	protected int max_it;
	protected double tol;
	protected double[] B, T, L, R, C; /* vectors of coeff matrix */
	protected int nx, ny;
	protected double dx, dy, da; /* mesh spacing parameters */

	protected double[][] phi, rho, b;
	final double EPS0 = 8.85418782e-12;

	public Solver(Field phi, Field rho, double dx, double dy) {
		/* convergence parameters */
		max_it = 5000;
		tol = 1e-4;

		/* set up dimensions */
		nx = phi.nx;
		ny = phi.ny;
		this.dx = dx;
		this.dy = dy;

		/* save fields */
		this.phi = phi.getData();
		this.rho = rho.getData();

		b = new Field(nx, ny).getData();

		/* allocate coefficient matrix: Bottom, Top, Left, Right, Center */
		B = new double[nx * ny];
		T = new double[nx * ny];
		L = new double[nx * ny];
		R = new double[nx * ny];
		C = new double[nx * ny];

		/* initialize coefficients */
		InitCoefficients();
	}

	/*
	 * InitBoundaries
	 */
	public void InitCoefficients() {
		int i, j, u;
		double v;

		double iddx = 1 / (dx * dx);
		double iddy = 1 / (dy * dy);

		/* set standard differencing on internal nodes */
		for (j = 1; j < ny - 1; j++)
			for (i = 1; i < nx - 1; i++) {
				u = j * nx + i;

				B[u] = iddy;
				T[u] = iddy;
				L[u] = iddx;
				R[u] = iddx;
				C[u] = -2 * (iddx + iddy);
			}

		/* set neumann boundaries along the edges */
		for (i = 0; i < nx; i++) {
			/* bottom */
			T[i] = 1 / dx;
			C[i] = -1 / dx;

			/* top */
			B[(ny - 1) * nx + i] = 1 / dx;
			C[(ny - 1) * nx + i] = -1 / dx;
		}

		for (j = 1; j < ny-1; j++) {
			/* left */
			R[j * nx] = 1 / dy;
			C[j * nx] = -1 / dy;

			/* right */
			L[j * nx + nx - 1] = 1 / dy;
			C[j * nx + nx - 1] = -1 / dy;
		}

		/* dirichlet boundary on top edge */
		for (i = nx / 3; i < 2 * nx / 3; i++) {
			u = (ny - 1) * nx + i;
			C[u] = 1;
			L[u] = R[u] = B[u] = T[u] = 0;
			phi[i][ny - 1] = 7;
		}

		/* dirichlet boundary on bottom edge */
		for (i = nx / 3; i < 2 * nx / 3; i++) {
			u = i;
			C[u] = 1;
			L[u] = R[u] = B[u] = T[u] = 0;
			phi[i][0] = 0.5;
		}

	}

	/* solves potential using Gauss-Seidel method */
	protected double update() {
		double tau;
		double norm = tol;
		int i, j;
		int it = 1; /*
					 * need to start with 1 so that we don't compute reside
					 * first time
					 */
		int u;

		/* recalculate RHS */
		for (j = 0; j < ny; j++)
			for (i = 0; i < nx; i++) {
				if (C[j * nx + i] == 1)
					b[i][j] = phi[i][j];
				else if (i == 0)
					b[i][j] = -10;
				else if (i == nx - 1)
					b[i][j] = 0;
				else if (j == 0)
					b[i][j] = 0;
				else if (j == ny - 1)
					b[i][j] = 0;
				else
					b[i][j] = -rho[i][j] / EPS0;
			}

		/* SOLVE POTENTIAL */
		while (it < max_it) {
			/* go only until ny-2 since the upper boundary is dirichlet */
			for (j = 0; j < ny; j++) {
				for (i = 0; i < nx; i++) {
					u = j * nx + i;

					tau = 0;
					if (j > 0)
						tau += B[u] * phi[i][j - 1];

					if (j < ny - 1)
						tau += T[u] * phi[i][j + 1];

					if (i > 0)
						tau += L[u] * phi[i - 1][j];

					if (i < nx - 1)
						tau += R[u] * phi[i + 1][j];

					phi[i][j] = (b[i][j] - tau) / C[u];
				}
			}

			/* check convergence */
			if (it % 50 == 0) {
				norm = 0;
				for (j = 0; j < ny; j++)
					for (i = 0; i < nx; i++) {
						u = j * nx + i;
						tau = 0;
						if (j > 0)
							tau += B[u] * phi[i][j - 1];

						if (j < ny - 1)
							tau += T[u] * phi[i][j + 1];

						if (i > 0)
							tau += L[u] * phi[i - 1][j];

						if (i < nx - 1)
							tau += R[u] * phi[i + 1][j];

						double R = phi[i][j] - (b[i][j] - tau) / C[u];
						norm += R * R;
					}
				norm = Math.sqrt(norm) / (nx * ny);
				System.out.printf(" Solver iteration %d, norm = %.2g\n", it, norm);

				if (norm < tol)
					break;
			}
			it++;
		}

		if (it >= max_it)
			System.err.printf(
					" !! Failed to converge in %d iterations, norm = %g\n", it,
					norm);
		else
			System.out.printf(" OK Converged in %d iterations, norm = %g\n", it,
					norm);

		return norm;
	}
}
