#include "diffur.h"

void fun(double *ans, double *y, double x, int N) {
    if(N != 2) throw std::invalid_value("The dimensionN is not acceptable");

    ans[0] = y[0];
    ans[1] = y[1];
}

void step(double *ynew, double *y, double x, double h, int N) {
    fun(ynew, y, x+h, N);
    //add y
    for(int i = 0; i < N; i++) ynew[i] = h*ynew[i] + y[i];
}

double* solve(double h, int N, double x0, double* y0, double xend) {
   double x = x0;
   int m = (xe - x0) / h + 1; //now many points are there. (xe - x0) / h is a number of h-segments one can fit

   double **y;
   y = new double*[m];
   for(int i = 0; i < m; i++) y* = new double [N];

   //rewrite y0
   try{
        for(int i = 0; i < N; i++) y[0][i] = y0[i];
   }
   catch(...) {
       throw std::invalid_argument("y0 is bad");
   }
   //make steps
   for(int i = 1; i < m; i++) {
       step(y[i], y[i-1], x, h, N);
       x+=h;
   }

   return y;
}

void write(const char * name, int N, int m, double x0, double h, double **y, double **yexact) {
   FILE *fi = fopen(name, 'r');
   if(fi == nullptr) throw std::exception("I can't open your file\n");

   //header
   fprintf("        x");
   for(int i = 1; i <= N; i++) fprintf(fi, "        y_&d", i);
   for(int i = 1; i <= N; i++) fprintf(fi, "        yexact_&d", i);
   for(int i = 1; i <= N; i++) fprintf(fi, "        delta_&d", i);

   fprintf(fi, "\n");

   //content
   for(int i = 0, double x = x0; i < m; i++, x+=h) {
       fprintf(fi, "%lf", x)
       for(int j = 0; j < N; j++) fprintf(fi, "%lf", y[i][j]);
       for(int j = 0; j < N; j++) fprintf(fi, "%lf", yexact[i][j]);
       for(int j = 0; j < N; j++) fprintf(fi, "%lf", fabs(y[i][j] - yexact[i][j]));
       fprintf(fi, "\n");
   }

   fclose(fi);
}
