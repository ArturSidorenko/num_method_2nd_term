#include "Header.h"

using namespace std;

int main() {
	double A;
	double x0 = 0;
	double y0 = 1;
	size_t m;
	double *y;
	double xend = 10;
	double h;

	FILE* fi;
	fopen_s(&fi, "output.txt", "w");
	if (fi == nullptr) throw std::bad_alloc();

	//header
	fprintf(fi, "Met.num  ");
	fprintf(fi, "    A         ");
	fprintf(fi, "    h         ");
	fprintf(fi, "     x        ");
	fprintf(fi, "     y        ");
	fprintf(fi, "  yexact        ");
	fprintf(fi, "delta       ");
	fprintf(fi, "delta *power  ");

	fprintf(fi, "\n");

	double norm;

	//convergence parameter
	for (int met = 1; met <= 5; met++) {
		for (A = 1; A <= 10; A ++) {
			for (h = 0.1; h > 10e-6; h *= 0.1) {
				solve(&y, &m, A, x0, y0, xend, h, met);
				//write down the last entity
				fprintf(fi, "%d    ", met);
				fprintf(fi, "% 12.4g  ", A);
				fprintf(fi, "% 12.4g  ", h);
				fprintf(fi, "% 12.4g  ", xend);
				fprintf(fi, "% 12.4g  ", y[m - 1]);
				fprintf(fi, "% 12.4g  ", exact(A, xend));
				norm = L2_norm(y, m, h, x0, A);
				fprintf(fi, "% 12.4g  ", norm);
				if(met == 3) fprintf(fi, "%9.2e ", norm / h / h);
				if(met != 3) fprintf(fi, "%9.2e ", norm / h );
				fprintf(fi, "\n");
				delete[] y;
			}
			fprintf(fi, "\n");
		}
		cout << "The method " << met << " has been processed\n";
		fprintf(fi, "\n\n");
	}

	fclose(fi);
	return 0;
}