#include "diffur.h"

using namespace std;

int main() {
    double **y;
    int N = 2;
    size_t m;
    double x0 = 0;
    double xend = 10;

    for(double h = 1; h > 1e-4; h *= 0.1) {

	double *y0 = new double[2];
	y0[0] = 1;
	y0[1] = 0;

	try {
                //write("output.txt", N, m, x0, h, y);
                cout << "The segment is [" << x0 << " , " << xend<<"], h = " << h << "\n";
                solve(&y, &m, h, N, x0, y0, xend);
                cout << "\t Old way.The Linf norm is equal to " << Linf_norm(N, m, x0, h, y) << "\n";
                for(size_t u = 0; u < m; u++) delete [] y[u];
                delete [] y;


                new_solve(&y, &m, h, N, x0, y0, xend);
                cout << "\t New way.The Linf norm is equal to " << Linf_norm(N, m, x0, h, y) << "\n";
                for(size_t u = 0; u < m; u++) delete [] y[u];
                delete [] y;


                newest_solve(&y, &m, h, N, x0, y0, xend);
                cout << "\t Newest way.The Linf norm is equal to " << Linf_norm(N, m, x0, h, y) << "\n";
	}
        catch (const std::invalid_argument &s) {
		std::cerr << "Don't worry. It is just some lousy error. " << s.what() << "\n";
	}
        catch (const std::exception &s) {
		std::cerr << "Don't worry. It is just some critical error. " << s.what() << "\n";
	}
	catch (...) {
		std::cerr << "Don't worry. It is just some unexpected error\n";
	}

        for(size_t u = 0; u < m; u++) delete [] y[u];
        delete [] y;
        delete [] y0;
    }

	return 0;
}
