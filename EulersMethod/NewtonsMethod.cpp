#include "Eulers.h"
#include "Gauss.h"
#include "NewtonsMethod.h"

vector<vector<double>> find_jacobian_1(const vector<double>& lambda)
{
    vector<vector<double>> jacobian_1(lambda.size(), vector<double>(lambda.size(), 0));
    jacobian_1.at(0).at(0) = diff_first_x1(lambda);
    jacobian_1.at(0).at(1) = diff_first_x2(lambda);
    jacobian_1.at(0).at(2) = diff_first_x3(lambda);
    jacobian_1.at(1).at(0) = diff_second_x1(lambda);
    jacobian_1.at(1).at(1) = diff_second_x2(lambda);
    jacobian_1.at(1).at(2) = diff_second_x3(lambda);
    jacobian_1.at(2).at(0) = diff_third_x1(lambda);
    jacobian_1.at(2).at(1) = diff_third_x2(lambda);
    jacobian_1.at(2).at(2) = diff_third_x3(lambda);
    return jacobian_1;
}

double find_d1(vector<double> x, vector<double> y, const vector<double>& lambdas)
{
    vector<double> elements = { abs(first_equasion(x, lambdas)),
                                     abs(second_equasion(x, lambdas)),
                                     abs(third_equasion(x, lambdas)) };

    return *max_element(elements.begin(), elements.end());
}

double find_d2(vector<double> x, vector<double> x1)
{
    size_t n = x.size();
    double max = 0;
    for (size_t i = 0; i < n; ++i) {
        if (abs(x1.at(i)) < 1) {
            if (abs(x1.at(i) - x.at(i)) > max) {
                max = abs(x1.at(i) - x.at(i));
            }
        }
        else {
            if (abs((x1.at(i) - x.at(i)) / x1.at(i)) > max) {
                max = abs((x1.at(i) - x.at(i)) / x1.at(i));
            }
        }
    }

    return max;
}

vector<vector<double>> fill_rectangle_matrix(const vector<vector<double>>& matrix, const vector<double>& const_vector) {
    vector<vector<double>> solvation_matrix;
    const size_t rows_num = matrix.size();
   
    for (size_t i = 0; i < rows_num; ++i) {
        vector<double> row;
        for (size_t j = 0; j <= rows_num; ++j) {
            if (j < rows_num) {
                row.push_back(matrix.at(i).at(j));
            }
            else {
                row.push_back(const_vector.at(i));
            }
        }
        solvation_matrix.push_back(row);
    }
    return solvation_matrix;
}



vector<double> solve_system_with_Newtons_method_three_equations(const double& e, vector<double>& x, const vector<double>& y, const vector<double>& lambdas, const double& h, const size_t& number_of_itterations )
{
    size_t n = x.size();
    size_t k_ = 1;
    double d1 = 1, d2 = 1;
    vector<double> F(3, 0);
    vector<vector<double>> J;
    vector<double> deltax;
    vector<double> x1(3, 0);
    while ((d1 > e || d2 > e) && k_ < number_of_itterations) {
        
        F.at(0) = -(first_equasion(x, lambdas));
        F.at(1) = -(second_equasion(x, lambdas));
        F.at(2) = -(third_equasion(x, lambdas));
     
        J = find_jacobian_1(lambdas);
        vector<vector<double>> rectangle = fill_rectangle_matrix(J, F);
        deltax = solve_system_with_gauss_method(rectangle, 0);
        for (size_t i = 0; i < n; ++i) {
            x1.at(i) = x.at(i) + deltax.at(i);
        }
        d1 = find_d1(x, y, lambdas);
        d2 = find_d2(x, x1);
        
        k_++;
        x = x1;
        if (k_ >= number_of_itterations) {
            cout << "ERROR: IER = 2" << endl;
        }
    }
    return x;
}
