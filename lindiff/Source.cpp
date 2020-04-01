#include "Header.h"

double exact(double A, double x) {
	return exp(-A * x);
}

void solve(double ** ans, size_t * size, double A, double x0, double y0, double xend, double h, int nummet)
{
	//make an array
	size_t m = (xend - x0) / h + 1;
	double *y;
	y = new double[m];

	y[0] = y0;

	double x = x0;
	
	if (nummet >= MET4) {
		y[1] = y0 - A * h;
		x += h;
		x += h;
		for (size_t i = 2; i < m; i++, x+=h) {
			y[i] = other_step(A, y[i - 1], y[i - 2], x, h, nummet);
		}
	}
	else {
		x += h;
		for (size_t i = 1; i < m; i++, x += h) {
			y[i] = step(A, y[i - 1], x, h, nummet);
		}
	}

	*ans = y;
	*size = m;
}

double step(double A, double y, double x, double h, int nummet)
{
	double ans;
	switch (nummet)
	{
	case MET1:
		ans = y - A * h * y;
		break;
	case MET2:
		ans = y / (1.0 + A * h);
		break;
	case MET3:
		ans = y * (1.0 - 0.5 * A * h) / (1.0 + 0.5 * A * h);
		break;
	default:
		throw std::invalid_argument("nummet");
		break;
	}
	return ans;
}

double other_step(double A, double y, double yold, double x, double h, int nummet)
{
	double ans;
	switch (nummet)
	{
	case MET4:
		ans = yold - 2 * A * h * y;
		break;
	case MET5:
		ans = 2.0 / 3.0 * (2 * y - 0.5 * yold - A * h * yold);
		break;
	default:
		throw std::invalid_argument("nummet");
		break;
	}
	return ans;
}

double L2_norm(double *y, size_t m, double h, double x0, double A) {
	double ans = 0;
	double x = x0;
	for (size_t i = 0; i < m; i++, x+=h) {
		ans += (y[i] - exact(A, x)) * (y[i] - exact(A, x)) * h;
	}
	return sqrt(ans);
}