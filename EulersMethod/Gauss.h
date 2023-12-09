#pragma once
#include "Eulers.h"
#include "Gauss.h"
#include "NewtonsMethod.h"

void swap_rows_for_not_null_diagonal(vector<vector<double>>& matrix);
vector<double> solve_system_with_gauss_method(vector<vector<double>>& solvation, const int& epsilon );