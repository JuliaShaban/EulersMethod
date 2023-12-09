#include "Eulers.h"
#include "Gauss.h"
#include "NewtonsMethod.h"


void swap_rows_for_not_null_diagonal(vector<vector<double>>& matrix) {
    double matrix_size = matrix.size();
    for (int i = 0; i < matrix_size; ++i) {
        if (abs(matrix.at(i).at(i)) <= round_error) {
            if (i > 0) {
                swap(matrix.at(i), matrix.at(i - 1));
            }
            else {
                swap(matrix.at(i), matrix.at(i + 1));
            }
        }
    }
}


vector<double> solve_system_with_gauss_method(vector<vector<double>>& solvation, const int& epsilon ) {
    swap_rows_for_not_null_diagonal(solvation);
    const double heigth = solvation.size();
    const double width = solvation.at(0).size();
   
    for (int i = 0; i < heigth; ++i) {  
        double aii = solvation.at(i).at(i);
         for (int j = i; j < width; ++j) {
            solvation.at(i).at(j) /= aii;
        }
       
        for (int j = i + 1; j < heigth; ++j) {
            double aji = solvation.at(j).at(i);
            for (int k = i + 1; k < width; ++k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
            solvation.at(j).at(i) = 0;
        }
    }
   
    for (int i = width - 2; i > 0; --i) {
       
        for (int j = i - 1; j >= 0; --j) {
            double aji = solvation.at(j).at(i);
            for (int k = width - 1; k >= i; --k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
        }
    }
    vector<double> answer;
    if (epsilon) {
        for (int i = 0; i < heigth; ++i) {
            answer.push_back((round((solvation.at(i).at(width - 1)) * pow(10, epsilon))) / pow(10, epsilon));
        }
    }
    else {
        for (int i = 0; i < heigth; ++i) {
            answer.push_back(solvation.at(i).at(width - 1));
        }
    }
    return answer;
}