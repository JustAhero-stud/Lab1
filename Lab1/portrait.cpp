#include "portrait.h"

Portrait::Portrait(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz)
{
	numEdges = (nz + 1) * (nx * (ny + 1) + ny * (nx + 1)) + (nx + 1) * (ny + 1) * nz;
	connections.resize(numEdges);

}

void Portrait::build(Matrix& A)
{
	build_connections();
	A.N = numEdges;
	A.IA.resize(A.N + 1);
	A.IA[0] = A.IA[1] = 0;

	for (int i = 2; i <= A.N; i++)
	{
		int col = A.IA[i - 1];
		A.IA[i] = col + connections[i - 1].size();
	}

	A.JA.resize(A.IA[A.N]);

	for (int i = 1, k = 0; i < A.N; i++)
		for (int j : connections[i])
		{
			A.JA[k] = j;
			k++;
		}

	A.DI.resize(A.N);
	A.AL.resize(A.IA[A.N]);
}

void Portrait::build_connections() {
	for (int q = 0; q < nz; q++) {
		for (int s = 0; s < ny; s++) {
			for (int p = 0; p < nx; p++) {
				for (int i = 0; i < 12; i++)
				{
					int ind1 = getGlobalNumber(nx, ny, p, s, q, i);
					for (int j = 1; j < 12; j++)
					{
						int ind2 = getGlobalNumber(nx, ny, p, s, q, j);
						if (ind2 > ind1)
							connections[ind2].insert(ind1);
					}
				}
			}
		}
	}
}