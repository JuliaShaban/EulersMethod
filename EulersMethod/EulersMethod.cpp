#include "Eulers.h"
#include "Gauss.h"
#include "NewtonsMethod.h"

int main()
{
    vector<double> start_x = { 10, 22, 9 };
    vector<double> lambdas = { -1, -1, -1 };
    print_vector(solve_system_of_three_differential_equations_with_Eulers_explicit_method(1,2, pow(10, -3)));
    cout << "\n\n\n\n";
    print_vector(solve_system_of_three_differential_equations_with_Eulers_implicit_method(start_x, lambdas, 1, 0 , pow(10, -5)));
    return 0;
}
