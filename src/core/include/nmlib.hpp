#pragma once
#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

#include "logger.hpp"

template <class Func, class Tuple, std::size_t... Is>
void apply_elemwise(Func&& f, const Tuple& tp, std::index_sequence<Is...>) {
    (f(std::get<Is>(tp)), ...);
}