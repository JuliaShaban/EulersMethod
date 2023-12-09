#include "Eulers.h"
#include "Gauss.h"
#include "NewtonsMethod.h"

double first_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((k - a) * x.at(1) * x.at(2)) / a;
}

double second_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((k + a) * x.at(0) * x.at(2)) / k;
}

double third_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((a - k) * x.at(0) * x.at(1)) / a;
}

double first_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) + 4 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(1) + (2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(2)) / 6 + (4 * lambda.at(0) + 2 * lambda.at(1)) / 6;
}

double second_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) * x.at(1) + (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) * x.at(2)) / 6 + (4 * lambda.at(0) - lambda.at(1) - 9 * lambda.at(2)) / 6;
}

double third_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) * x.at(1) + (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) * x.at(2)) / 6 + (4 * lambda.at(0) - lambda.at(1) + 9 * lambda.at(2)) / 6;
}

double diff_first_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) + 4 * lambda.at(1)) / 6;
}

double diff_first_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_first_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_second_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_second_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) / 6;
}

double diff_second_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) / 6;
}

double diff_third_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_third_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) / 6;
}

double diff_third_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) / 6;
}
double find_step(double e, double hmax, vector<double>& answer) {
    double h1, h2, h3, hmin;
    h1 = e / (abs(answer.at(0)) + (e / hmax));
    h2 = e / (abs(answer.at(1)) + (e / hmax));
    h3 = e / (abs(answer.at(2)) + (e / hmax));
    if (h1 < h2) {
        hmin = h1;
    }
    else {
        hmin = h2;
    }
    if (h3 < hmin) {
        hmin = h3;
    }
    return hmin;
}

vector<double> solve_system_of_three_differential_equations_with_Eulers_explicit_method(const int& a, const int& k , const double& e ) {
    double t = 0, T = 1;
    vector<double> answer(3, 1);
    double hmax = 1, h;
    vector<double> answer_1(3, 0);
    size_t i = 0;

    do {
        i++;
        answer_1.at(0) = first_equasion_1(answer, a, k);
        answer_1.at(1) = second_equasion_1(answer, a, k);
        answer_1.at(2) = third_equasion_1(answer, a, k);
        h = find_step(e, hmax, answer_1);
        answer.at(0) += h * answer_1.at(0);
        answer.at(1) += h * answer_1.at(1);
        answer.at(2) += h * answer_1.at(2);
        t += h;
        cout <<"t = "<< setprecision(6)  << t << ";" << setw(10) << answer.at(0) << ";" << setw(5) << answer.at(1) << ";" << setw(5) << answer.at(2) << endl;
    } while (t < T);
    return answer;
}

double find_step_1(double h, double e, const vector<double>& ek, bool v)
{
    vector<double> hh(3, 0);
    if (v) {
        for (size_t i = 0; i < 3; ++i) {
            hh.at(i) = pow(e / abs(ek.at(i)), 0.5) * h;
        }
    }
    else {
        for (size_t i = 0; i < 3; ++i) {
            if (abs(ek.at(i)) > e) {
                hh.at(i) = h / 2;
            }
            if (e / 4 < abs(ek.at(i)) && abs(ek.at(i)) <= e) {
                hh.at(i) = h;
            }
            if (abs(ek.at(i)) < e / 4) {
                hh.at(i) = 2 * h;
            }
        }
    }
    return *min_element(hh.begin(), hh.end());
}


vector<double> find_new_eps(vector<double> uk_1, vector<double> u, vector<double> uk1, double h, double hk_1)
{
    vector<double> ek(3, 0);
    ek.at(0) = -(h / (h + hk_1)) * (uk1.at(0) - u.at(0) - (h * (u.at(0) - uk_1.at(0))) / (hk_1));
    ek.at(1) = -(h / (h + hk_1)) * (uk1.at(1) - u.at(1) - (h * (u.at(1) - uk_1.at(1))) / (hk_1));
    ek.at(2) = -(h / (h + hk_1)) * (uk1.at(2) - u.at(2) - (h * (u.at(2) - uk_1.at(2))) / (hk_1));
    return ek;
}

vector<double> solve_system_of_three_differential_equations_with_Eulers_implicit_method(const vector<double>& start_x, const vector<double>& lambdas, const double& T, const bool& mode , const double& e ) {
    double t = 0;
    vector<double> answer = start_x;
    double hmax = 1;
    double hmin = 0.01;
    double h = hmin;
    double hk1 = hmin;
    double hk_1 = hmin;
    double tk;
    vector<double> uu(3, 0);
    vector<double> uk1(3, 0);
    vector<double> uk_1(3, 0);

    vector<double> ek(3, 0);
    int i = 0, z;
    do {
        do {
            z = 0;
            i++;
            tk = t + h;
            uk1 = solve_system_with_Newtons_method_three_equations(e, uk1, answer, lambdas, h, 500);
            ek = find_new_eps(uk_1, answer, uk1, h, hk_1);
            if (abs(ek.at(0)) > e || abs(ek.at(1)) > e || abs(ek.at(2)) > e) {
                h /= 2;
                hk1 = h;
                uk1.at(0) = answer.at(0);
                uk1.at(1) = answer.at(1);
                uk1.at(2) = answer.at(2);
                z = 1;
            }
        } while (z == 1);
  
        hk1 = find_step_1(h, e, ek, mode);
        if (hk1 > hmax) {
            hk1 = hmax;
        }
        uk_1.at(0) = answer.at(0);
        uk_1.at(1) = answer.at(1);
        uk_1.at(2) = answer.at(2);
        answer.at(0) = uk1.at(0);
        answer.at(1) = uk1.at(1);
        answer.at(2) = uk1.at(2);
        hk_1 = h;
        h = hk1;
        t = tk;
        cout << setprecision(6) << "t = " << fixed  << t << ";" << setw(15) << answer.at(0) << ";" << setw(5) << answer.at(1) << ";" << setw(5) << answer.at(2) << endl;
    } while (t < T);
    return answer;
}

void print_matrix(const vector<vector<double>>& matrix) {
    const double rows_num = matrix.size();
    const double column_num = matrix.at(0).size();
    for (int i = 0; i < rows_num; ++i) {
        for (int j = 0; j < column_num; ++j) {
            cout << setw(outwidth) << matrix.at(i).at(j);
        }
        cout << endl;
    }
}


void print_vector(const vector<double>& const_vector) {
    double length = const_vector.size();
    for (int i = 0; i < length; ++i) {
        cout << setw(outwidth) << const_vector.at(i);
    }
    cout << endl;
}

