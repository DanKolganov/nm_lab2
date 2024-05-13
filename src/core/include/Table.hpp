#pragma once

#include <cstddef>
#include <vector>
struct TableRow {
    double x;
    double v;
    double u; // real solv
    double eps; // epsilon in point ||ui-vi||
};

class Table {
using rows = std::vector<TableRow>;
    rows _data;
public:
    Table() = default;
    virtual ~Table() = default;

    Table(const size_t& n) {
        rows temp;
        temp.reserve(n);
        _data = temp;
    }

    Table(rows&& data) : _data(data) {}

    virtual void reserve(const size_t& n) {
        _data.reserve(n);
    }

    virtual void append(const TableRow& row) {
        _data.push_back(row);
    };

    virtual TableRow at(const size_t& indx) { 
        return _data.at(indx);
    }

    TableRow& operator[](const size_t& indx) { 
        return _data[indx];
    }

    virtual size_t size() const noexcept {
        return _data.size();
    }
};