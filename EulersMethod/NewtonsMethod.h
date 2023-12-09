#pragma once
vector<vector<double>> find_jacobian_1(const vector<double>& lambda);
double find_d1(vector<double> x, vector<double> y, const vector<double>& lambdas);
double find_d2(vector<double> x, vector<double> x1);
vector<double> solve_system_with_Newtons_method_three_equations(const double& e, vector<double>& x, const vector<double>& y, const vector<double>& lambdas, const double& h, const size_t& number_of_itterations );
vector<vector<double>> fill_rectangle_matrix(const vector<vector<double>>& matrix, const vector<double>& const_vector);