#include "eigen.h"

const double PI = 4 * atan(1);

double eigenvalue(unsigned N, unsigned q) {

    double h = 2./(2* N - 1);
    return 2. / h / h * (1 - cos( (q + 0.5) / (N - 0.5) * PI ) );
}

double eigenvector(unsigned N, unsigned q, unsigned k) {
    double h = 2. / (2 * N - 1);

    return sqrt(2) * sin((q + 0.5) * PI * h * k);
}

double scal_product(unsigned N, unsigned p, unsigned q) {
    double sum = 0, h = 2. /(2 * N - 1);
    for(unsigned i = 1; i < N; i++) {
        sum += eigenvector(N, p, i) * eigenvector(N, q, i) * h;
    }
    return sum;
}

double left_hand(unsigned N, unsigned q, unsigned k) {
    double h = 2. /(2 * N - 1);

    double numerator, denominator;

    numerator = eigenvector(N, q, k+1) - 2 * eigenvector(N, q, k) + eigenvector(N, q, k-1);
    denominator = h * h;

    return - numerator / denominator;
}

double right_hand(unsigned N, unsigned q, unsigned k) {
    return eigenvalue(N, q) * eigenvector(N, q, k);
}

double eq_res(unsigned N, unsigned q) {
    double res = 0, cur = 0;
    for(unsigned j = 1; j < N; j++) {
        cur = fabs((left_hand(N, q, j) - right_hand(N, q, j)) / eigenvalue(N, q));
        res = std::max(res, cur);
    }

    //include boundary conditions
    cur = fabs(eigenvector(N, q, 0) / eigenvalue(N, q));
    res = std::max(res, cur);
    cur = fabs((eigenvector(N, q, N-1) - eigenvector(N, q, N)) / eigenvalue(N, q));
    res = std::max(res, cur);

    return res;
}

double eq_res_all(unsigned N) {
    double res = 0, cur = 0;
    for(unsigned q = 0; q <= (N - 2); q++) {
        cur = eq_res(N, q);
        res = std::max(cur, res);
    }

    return res;
}

double orthogonality_check(unsigned N) {
    double res = 0, cur = 0;
    for(unsigned q = 0; q <= (N - 2); q++) {
        for(unsigned p = 0; p <=(N-2); p++) {
            cur = scal_product(N, p, q);
            if(p == q) cur-=1.;
            res = std::max(res, fabs(cur));
        }
    }
    return res;
}

void print_eigenv(unsigned N) {
    FILE *f;
    f = fopen("eigenvals.txt", "w");


    fprintf(f, "Number      Discrete        Continious\n");

    for(unsigned i = 0; i <= (N-2); i++) {
        fprintf(f, "%d       ", i);
        fprintf(f, "%9.5lf   ", eigenvalue(N, i));
        fprintf(f, "%9.5lf   ", (PI * 0.5 + PI*(i)) * (PI * 0.5 + PI*(i)));
        fprintf(f, "\n");
    }

    fclose(f);
}

void print_gramian(unsigned N) {
    FILE *f;
    f = fopen("Gramian.txt", "w");


    for(unsigned i = 0; i <= (N-2); i++) {
        for(unsigned j = 0; j <= (N-2); j++) {
            fprintf(f, " %9.5e ", scal_product(N, i, j));
        }
        fprintf(f, "\n");
    }

    fclose(f);

}
void print_vector(unsigned N, unsigned q) {
    FILE *f;
    f = fopen("Vec_disc.txt", "w");
    const double h = 2./(2*N-1);

    for(unsigned j = 0; j <=N; j++) {
        fprintf(f, " %9.5lf  ", h*j);
        fprintf(f, " %9.5lf  ", eigenvector(N, q, j));
        fprintf(f, " %9.5lf  ", sqrt(2) * sin(j*h * (PI*0.5 + PI*q)));
        fprintf(f, "\n");
    }

    fclose(f);
    FILE *g;
    g = fopen("Vec_cont.txt", "w");
    for(unsigned j = 0; j <=300; j++) {
        fprintf(g, " %9.5lf  ", j * (1./300));
        fprintf(g, " %9.5lf  ", sqrt(2) * sin(j*(1./300) * (PI*0.5 + PI*q)));
        fprintf(g, "\n");
    }
    fclose(g);

}
