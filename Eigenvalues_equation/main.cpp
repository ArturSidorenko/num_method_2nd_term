#include "eigen.h"

using namespace std;

int main() {
    unsigned N;

    cout << "Enter N: ";
    cin >> N;
    cout<<"\n";

    cout << "The maximal relative residual of equation (boundary conditions are included): " << eq_res_all(N) <<"\n\n";

    cout << "The difference between Gramian and identity matrices: " << orthogonality_check(N) << "\n\n";

    return 0;

}
