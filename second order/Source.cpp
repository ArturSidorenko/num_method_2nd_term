#include "Header.h"

extern double PI = 4 * atan(1);

task::task(unsigned N, double b, myfunc right_hand): N_(N)
{
	size_ = N - 1;
	h_ =  2.0 /(2 * N - 1);
	isconst_ = true;
	f_.resize(size_);
	bconst_ = b;
	for (size_t i = 1; i < N; i++) f_[i-1] = right_hand( h_ * i );
}

task::task(unsigned N, myfunc b, myfunc right_hand): N_(N)
{
	isconst_ = false;
	size_ = N - 1;
	h_ = 2.0 / (2 * N - 1);
	b_.resize(size_);
	f_.resize(size_);
	for (size_t i = 1; i < N; i++) {
		b_[i - 1] = b( h_ * i ) ;
	}
	for (size_t i = 1; i < N; i++) f_[i-1] = right_hand( h_ * i ) ;
}

double task::getb() const
{
	if (!isconst_) throw std::invalid_argument("b is not constant");
	return bconst_;
}

double task::getb(size_t i) const
{
	if (isconst_) return bconst_;
	return b_[i];
}

lin_sys::lin_sys(size_t m): m_(m)
{
	a_.resize(m + 2, 0);
	b_.resize(m + 2, 0);
	c_.resize(m + 2, 0);
	f_.resize(m + 2, 0);
}

double basis(unsigned N, unsigned k, unsigned p)
{
	double h = 2.0 / (2 * N - 1);
	return sqrt(2) * sin(PI * (k - 0.5) * p * h);
}

double eigenval(unsigned N, unsigned k)
{
	double h = 2.0 / (2 * N - 1);
	return 2.0 / h /h * (1 - cos(PI * (k - 0.5) * h));
}


void progonka(myvec & v, const lin_sys &s)
{
	size_t m = s.m_;
	const myvec &a = s.a_, &b = s.b_, &c = s.c_, &f = s.f_;

	v.clear();
	v.resize(m + 1, 0);
	myvec alpha(m + 2), beta(m + 2);

	//direct step
	alpha[0] = beta[0] = 0;
	for (size_t i = 0; i <= m; i++) {
		alpha[i + 1] = (-b[i]) / (c[i] - (-a[i]) * alpha[i]);
		beta[i + 1] = (f[i] + (-a[i]) * beta[i]) / (c[i] - (-a[i]) * alpha[i]);
	}

	//reverse step
	v[m] = beta[m + 1];
	for (size_t i = m; i > 0; i--) v[i - 1] = v[i] * alpha[i] + beta[i];
}

void prepare_progonka(myvec &v, const task &t) {
	const double h = t.h_;
	const double w = 1.0 / h / h;
	
	size_t m = t.size_ - 1;
	lin_sys s(m);

	for (size_t i = 0; i < m; i++) s.c_[i] = 2.0 * w  + t.getb(i);
	s.c_[m] = w + t.getb(m);
	for (size_t i = 1; i <= m; i++) s.a_[i] = -w ;
	for (size_t i = 0; i <  m; i++) s.b_[i] = -w ;
	for (size_t i = 0; i <= m; i++) s.f_[i] = t.f_[i];


	progonka(v, s);
}


void fourier(myvec &v, const task &t) {
	if (!t.isconst_) throw std::invalid_argument("b is not a constant");
	double b = t.bconst_;
	size_t m = t.size_;
	myvec c(m);
	v.clear();
	v.resize(m);

	for (size_t i = 1; i <= m; i++) {
		//scalar product
		c[i - 1] = 0;
		for (size_t p = 1; p <= m; p++) {
			c[i - 1] += t.f_[p-1] * basis(t.N_, i, p) * t.h_;
		}
		c[i - 1] /= eigenval(t.N_, i) + b;
	}

	//now we calculate the answer
	for (size_t p = 1; p <= m; p++) {
		v[p - 1] = 0;
		for (size_t i = 1; i <= m; i++) v[p - 1] += c[i - 1] * basis(t.N_, i, p);
	}
}

double L2_norm(const myvec & a, const myvec & b, size_t m, double h)
{
	double ans = 0;
	for (size_t i = 0; i < m; i++) ans += (b[i] - a[i]) * (b[i] - a[i]) * h;
	return sqrt(ans);
}