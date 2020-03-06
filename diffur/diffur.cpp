#include "diffur.h"

void fun(double *ans, const double *y, double x, int N) {
    if(N != 2) throw std::invalid_argument("The dimensionN is not acceptable");

    ans[0] = y[1] + x - x;
    ans[1] = -y[0];
}

void fexact(double *y, double x, int N) {
	if (N != 2) throw std::invalid_argument("fexact");


	y[0] =  cos(x);
	y[1] = -sin(x);

}

void step(double *ynew, const double *y, double x, double h, int N) {
    fun(ynew, y, x, N);
    //add y
    for(int i = 0; i < N; i++) ynew[i] = h * ynew[i] + y[i];
}

void new_step(double *ynew, const double *y, const double *yold, double x, double h, int N) {
    fun(ynew, y, x, N);

    for(int i = 0; i < N; i++) ynew[i] = 2 * h * ynew[i] + yold[i];
}

void solve(double ***ans, size_t *points, double h, int N, double x0, const double* y0, double xend) {
   double x = x0;
   size_t m = (xend - x0) / h + 1; //now many points are there. (xe - x0) / h is a number of h-segments one can fit


   auto y = new double*[m];
   for(size_t i = 0; i < m; i++) y[i] = new double [N];

   //rewrite y0
   try{
        for(int i = 0; i < N; i++) y[0][i] = y0[i];
   }
   catch(...) {
       throw std::invalid_argument("the y0 is bad");
   }
   //make steps
   for(size_t i = 1; i < m; i++) {
       step(y[i], y[i-1], x, h, N);
       x+=h;

   }

   *ans = y;
   *points = m;

}

void new_solve(double ***ans, size_t *points, double h, int N, double x0, const double* y0, double xend) {
   double x = x0;
   size_t m = (xend - x0) / h + 1; //now many points are there. (xe - x0) / h is a number of h-segments one can fit


   auto y = new double*[m];
   for(size_t i = 0; i < m; i++) y[i] = new double [N];

   //rewrite y0
   try{
        for(int i = 0; i < N; i++) y[0][i] = y0[i];
   }
   catch(...) {
       throw std::invalid_argument("the y0 is bad");
   }

   //calculate y1
   fexact(y[1], x + h, N);

   //make steps
   for(size_t i = 2; i < m; i++) {
       new_step(y[i], y[i-1], y[i-2], x, h, N);
       x+=h;
   }

   *ans = y;
   *points = m;
}



void write(const char *name, int N, int m, double x0, double h, double **y) {
   FILE* fi;
   fi = fopen(name, "w");
   if(fi == nullptr) throw std::bad_alloc();

   //header
   fprintf(fi, "      x  ");
   for(int i = 1; i <= N; i++) fprintf(fi, "     y_%d      ", i);
   for(int i = 1; i <= N; i++) fprintf(fi, "yexact_%d      ", i);
   for(int i = 1
; i <= N; i++) fprintf(fi, "delta_%d       ", i);

   fprintf(fi, "\n");

   //content
   double x = x0;
   double *yexact = new double[N];
   for(int i = 0; i < m; i++, x+=h) {
	   fprintf(fi, "%9.5lf    ", x);
	   fexact(yexact, x, N);
       for(int j = 0; j < N; j++) fprintf(fi, "%9.5lf    ", y[i][j]);
	   for (int j = 0; j < N; j++) fprintf(fi, "%9.5lf    ", yexact[j]);
       for(int j = 0; j < N; j++) fprintf(fi, "%9.5lf    ", fabs(y[i][j] - yexact[j]));
       fprintf(fi, "\n");
	   
   }

   delete[] yexact;
   fclose(fi);
}

double Linf_norm(int N, int m, double x0, double h, double **y) {
    double ans = 0;
    double *yexact = new double[N];
    double x = x0;
    for(int i = 0; i < m; i++, x+=h) {
        for(int j = 0; j < N; j++) {
            fexact(yexact, x, N);
            ans = std::max(ans, fabs(y[i][j] - yexact[j]));
        }
    }
    delete [] yexact;
    return ans;
}

