//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_PREPROCESS_TIME_WINDOWS_H
#define EJOR2019_PREPROCESS_TIME_WINDOWS_H

#include <goc/goc.h>

namespace ejor2019
{
// Takes a JSON instance of the vehicle routing problems with the following attributes:
//	- digraph
//	- travel_times
//	- time_windows
//	- start_depot
//	- end_depot
// Assumes preprocess_travel_times was called (i.e. instance has no service nor waiting times).
// Shrinks time windows [a_i, b_i] according to the techinques introduced in:
//	Desrosiers, J., Dumas, Y., Solomon, M. M., & Soumis, F. (1995).
// and removes infeasible arcs.
// 	* is_profitable: indicates is the instance is a profitable problem. In that case we need to remove some
// 					 preprocessings.
void preprocess_time_windows(nlohmann::json& instance, bool is_profitable);
} // namespace ejor2019

#endif //EJOR2019_PREPROCESS_TIME_WINDOWS_H
