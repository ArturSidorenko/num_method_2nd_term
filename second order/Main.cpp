#include "Header.h"

using namespace std;
extern double PIPI = 4 * atan(1);

double sol1(double x) {
	return x * exp(x) - 2 * exp(1) * x ;
}

double f1(double x) {
	double a, c, b = 0;
	a = 2 * exp(x) + x * exp(x);
	c = b * sol1(x);
	return c - a;
}

double f2(double x) {
	double a, c, b = 10000;
	a = 2 * exp(x) + x * exp(x);
	c = b * sol1(x);
	return c - a;
}

double myb(double x) {
	//return sqrt(x + 0.42) * sin(x + 0.4242) * exp(42 * x);
	return 100;
}

double f3(double x) {
	double a, c, b = myb(x);
	a = 2 * exp(x) + x * exp(x);
	c = b * sol1(x);
	return c - a;
}

void test(const char *name, myfunc f, double b) {

	FILE* fi;
	fopen_s(&fi, name, "w");
	if (fi == nullptr) throw std::bad_alloc();


	myvec ans, ans2, exact;
	fprintf(fi, "b = %g\n\n", b);
	fprintf(fi, "Number of segments   ");
	fprintf(fi, " L2 norm between Exact and Progonka   ");
	fprintf(fi, " L2 norm between Exact and Fourier    ");
	fprintf(fi, " L2 norm between Progonka and Fourier  ");
	fprintf(fi, " Factor(L2 norm / h ** 2)\n\n");
	for (size_t N = 5; N < 3000; N *= 2) {
		task t(N, b, f);

		exact.resize(N - 1);
		for (size_t i = 1; i < N; i++) exact[i - 1] = sol1(t.h_ * i);

		//solve via progonka
		prepare_progonka(ans, t);

		fprintf(fi, "%8d                 ", N);
		fprintf(fi, "%12.4g                        ", L2_norm(exact, ans, t.size_, t.h_));

		fourier(ans2, t);

		fprintf(fi, "%12.4g                           ", L2_norm(exact, ans2, t.size_, t.h_));
		fprintf(fi, "%12.4g                  ", L2_norm(ans2, ans, t.size_, t.h_));
	
		fprintf(fi, "% 12.4g \n", L2_norm(exact, ans, t.size_, t.h_) / t.h_ / t.h_);
	}
}

void test(const char *name, myfunc f, myfunc b) {

	FILE* fi;
	fopen_s(&fi, name, "w");
	if (fi == nullptr) throw std::bad_alloc();


	myvec ans, ans2, exact;
	fprintf(fi, "b is not a constant.\n\n");
	fprintf(fi, "Number of segments   ");
	fprintf(fi, " L2 norm between Exact and Progonka   ");
	fprintf(fi, " Factor(L2 norm / h ** 2)\n\n");
	for (size_t N = 5; N < 3000; N *= 2) {
		task t(N, b, f);

		exact.resize(N - 1);
		for (size_t i = 1; i < N; i++) exact[i - 1] = sol1(t.h_ * i);

		//solve via progonka
		prepare_progonka(ans, t);

		fprintf(fi, "%8d                 ", N);
		fprintf(fi, "%12.4g                        ", L2_norm(exact, ans, t.size_, t.h_));


		fprintf(fi, "% 12.4g \n", L2_norm(exact, ans, t.size_, t.h_) / t.h_ / t.h_);
	}
}

int main() {

	test("output_zero.txt", &f1, 0.0);
	test("output_10k.txt", &f2, 10000);
	test("output_var.txt", &f3, myb);

	system("pause");
	return 0;
}