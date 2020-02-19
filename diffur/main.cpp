#include "diffur.h"

using namespace std;

int main() {
    double **y;
    int N = 2;
    double x0 = 0;
    double xend = 1;
    double h = 0.001;

	double *y0 = new double[2];
	y0[0] = 1;
	y0[1] = 0;

	try {

		y = solve(h, N, x0, y0, xend);

		write("output.txt", N, x0, xend, h, y);
	}
	catch (std::invalid_argument s) {
		std::cerr << "Don't worry. It is just some lousy error. " << s.what() << "\n";
	}
	catch (std::exception s) {
		std::cerr << "Don't worry. It is just some critical error. " << s.what() << "\n";
	}
	catch (...) {
		std::cerr << "Don't worry. It is just some unexpected error\n";
	}

	return 0;
}
