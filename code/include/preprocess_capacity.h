//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_PREPROCESS_CAPACITY_H
#define EJOR2019_PREPROCESS_CAPACITY_H

#include <goc/goc.h>

namespace ejor2019
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- capacity
//	- demands
// Removes arcs (i, j) such that q_i+q_j > Q.
void preprocess_capacity(nlohmann::json& instance);
} // namespace ejor2019

#endif //EJOR2019_PREPROCESS_CAPACITY_H
