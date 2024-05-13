#include "diffschemes.hpp"

#include "Table.hpp"
#include "nmlib.hpp"
#include "nmlinalg.hpp"
#include "logger.hpp"
#include <cmath>
#include <cstddef>
#include <math.h>

NMbalance::~NMbalance() {}

void NMbalance::eval() {
    const double& x = _net->first;  // x coord of node in net
    const double& step = _net->step;
    const double& n = _net->n;
    const double stepsq = step * step;

    LOG_DEBUG_CLI("net: ", _net->n, _net->step, _net->first, _net->last);
    // evaluating the coefficients for linear system (tridiagonal matrix algorithm(TMA) aka "progonka")
    auto cA = [&](const size_t& i) -> double {
        return _ca(x + (i - 1) * step, x + i * step, step) / stepsq;;
    };
    auto cB = [&](const size_t& i) -> double {
        return -((_ca(x + (i - 1) * step, x + i * step, step) + _ca(x + (i)*step, x + (i + 1) * step, step)) / stepsq +
                 _cd(x + (i - 0.5) * step, x + (i + 0.5) * step, step));
    };
    auto cC = [&](const size_t& i) -> double {
        return _ca(x + i * step, x + (i + 1) * step, step) / stepsq;
    };

    auto cF = [&](const size_t& i) -> double {
        return -_cphi(x + (i - 0.5) * step, x + (i + 0.5) * step, step);
    };

    LOG_DEBUG_CLI("Evalutating matrix coeffs...");
    matrix[0][0] = 0.l;
    matrix[0][1] = 1.l;
    matrix[0][2] = 0.l;
    rhs[0] = _mu1;

    for (size_t i = 1; i < n; i++) {
        matrix[i][0] = cA(i);
        matrix[i][1] = cB(i);
        matrix[i][2] = cC(i);
        rhs[i] = cF(i);
    }

    matrix[n][0] = 0.l;
    matrix[n][1] = 1.l;
    matrix[n][2] = 0.l;
    rhs[n] = _mu2;
    // coefficient evaluated and wrote in matrix and rhs vars
    // now we need to solve linear system w/ TMA

    LOG_DEBUG_CLI("Applying TMA algorithm...");
    TridiagonalAlg(matrix, vars, rhs, n + 1);

    for (size_t i = 0; i < vars.size(); i++) {
    _result->append({(x + step * i), vars[i],2,std::abs(vars[i]-2)});
    }

    // actually we need to put solution in some resultTable like object, but for some time it will be more... reasonable
}