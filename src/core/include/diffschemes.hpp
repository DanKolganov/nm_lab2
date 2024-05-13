#pragma once

#include "Table.hpp"
#include "logger.hpp"
#include "nmlib.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

struct Uniform1DNet {
    /* Uniform1DNet is uniform net with nodes: 
    xi = first + i * step for i=0,n and step = (last - first / n) */

    double first;
    double last;
    size_t n;

    double step;

    Uniform1DNet() : first(0.l), last(0.l), n(0), step(0.l) {}

    Uniform1DNet(const double& _a, const double& _b, const double& _n) {
        first = _a;
        last = _b;
        n = _n;

        step = (last - first) / n;
    }
};

class NMbalance {
    using func = std::function<double(double, double, double)>;

    Uniform1DNet* _net;
    Table* _result;

    /* coeff functions for diff scheme */
    func _ca;
    func _cd;
    func _cphi;
    double _mu1, _mu2;

    std::vector<std::vector<double>> matrix; /* matrix of linear equations system, with C-style layout
                                (probably tridiagonal so it can looks like this (v11, v12, v13, v21, v22 .. etc)) */
    std::vector<double> rhs;                 /* linear system rhs */

    std::vector<double> vars; /* linear system vars */

public:
    NMbalance() = default;
    ~NMbalance();

    NMbalance(Uniform1DNet* net, Table* tbl, func&& ca, func&& cd, func&& cphi, const double& mu1, const double& mu2)
        : _net(net), _result(tbl), _ca(ca), _cd(cd), _cphi(cphi), _mu1(mu1), _mu2(mu2) {
        rhs.resize(net->n + 1);
        matrix.resize((net->n + 1));
        for (auto& v : matrix) {
            v.resize(3);
        }
        vars.resize(net->n + 1);
    }

    void eval();

    Table* getTable() {
        return _result;
    }

    void setNet(Uniform1DNet* net) {
        _net = net;
    }
};