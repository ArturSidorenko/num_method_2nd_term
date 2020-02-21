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
