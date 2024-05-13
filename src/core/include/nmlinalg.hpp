#pragma once

#include <cstddef>
#include "logger.hpp"
template <class T, class U>
void TridiagonalAlg(T& coeffs, U& x, U& rhs, const size_t& n) {
    /* coeffs is an object, which contains a TMA coefficients 
       and can that coefficients can be accessed in form coeffs[i][j] where j =0,1,2 
       x is an object, which represent a vector of linear system variables which can 
       be accessed in form x[i] (x[0] = x0, x[1] = x1...)
       n is number of variables in linear system 
       rhs is an object, which represent a vector of linear system right-hand-sides which can 
       be accessed in form rhs[i] (rhs[0] = F1 (?)) */

    double temp;
    for (size_t i = 1; i < n; i++) { 
        temp = coeffs[i][0] / coeffs[i - 1][1];
        coeffs[i][1] = coeffs[i][1] - temp * coeffs[i - 1][2];
        rhs[i] = rhs[i] - temp * rhs[i - 1];
    }

    x[n - 1] = rhs[n - 1] / coeffs[n - 1][1];

    for (long long int i = n - 2; i >= 0; i--) {
        x[i] = (rhs[i] - coeffs[i][2] * x[i + 1]) / coeffs[i][1];
    }
}