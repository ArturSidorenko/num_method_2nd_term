#include "eigen.h"

using namespace std;

int main() {
    unsigned N;

    cout << "Enter N: ";
    cin >> N;
    cout<<"\n";

    cout << "The maximal relative residual of the equations (the boundary conditions are included): " << eq_res_all(N) <<"\n\n";

    cout << "The difference between the Gramian and the identity matrices: " << orthogonality_check(N) << "\n\n";

    print_eigenv(N);
    print_gramian(N);
    print_vector(N, 2);

    return 0;

}
