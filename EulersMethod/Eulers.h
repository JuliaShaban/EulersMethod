#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <algorithm>
using namespace std;

const size_t outwidth = 20;
const size_t precision = 10;
const double round_error = 0.00001;
double first_equasion_1(const vector<double>& x, const int& a, const int& k);
double second_equasion_1(const vector<double>& x, const int& a, const int& k);
double third_equasion_1(const vector<double>& x, const int& a, const int& k);
double first_equasion(const vector<double>& x, const vector<double>& lambda);
double second_equasion(const vector<double>& x, const vector<double>& lambda);
double third_equasion(const vector<double>& x, const vector<double>& lambda);
double diff_first_x1(const vector<double>& lambda);
double diff_first_x2(const vector<double>& lambda);
double diff_first_x3(const vector<double>& lambda);
double diff_second_x1(const vector<double>& lambda);
double diff_second_x2(const vector<double>& lambda); 
double diff_second_x3(const vector<double>& lambda);
double diff_third_x1(const vector<double>& lambda);
double diff_third_x2(const vector<double>& lambda);
double diff_third_x3(const vector<double>& lambda);
double find_step(double e, double hmax, vector<double>& answer);
vector<double> solve_system_of_three_differential_equations_with_Eulers_explicit_method(const int& a , const int& k, const double& e );
double find_step_1(double h, double e, const vector<double>& ek, bool v);
vector<double> find_new_eps(vector<double> uk_1, vector<double> u, vector<double> uk1, double h, double hk_1);
vector<double> solve_system_of_three_differential_equations_with_Eulers_implicit_method(const vector<double>& start_x, const vector<double>& lambdas, const double& T, const bool& mode , const double& e );
void print_matrix(const vector<vector<double>>& matrix);
void print_vector(const vector<double>& const_vector);