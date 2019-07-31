//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_TDFI_H
#define EJOR2019_TDFI_H

#include <goc/goc.h>
#include "vrp_instance.h"

namespace ejor2019
{
// The TDFI are the following:
// Let k \in V - {o, d} be a vertex, F \subseteq {(i,k,m) \in \hat{A}}
// 		\sum_{(i,k,m) \in F} x_ikm <= \sum_{(k,j,m) \in \delta^+(F)} x_kjm
//
// The separation routine does the following:
//	Let AI(e) = "arrival interval of arc e", DI(e) = "departure interval of arc e".
// 	1) Fixes a vertex k \in V - {o,d}.
// 	2) Sets I = {e_1, ..., e_r}, O = {f_1, ..., f_s} the sets of arcs entering k and leaving k respectively (with z*[e_i], z*[f_i] > 0).
//	3) Sorts I by the start of their arrival interval, and sorts O by the start of their departure interval.
// 	4) Precompute next[i] = min { j: min(DI(f_j)) > max(AI(e_i)) or j = s + 1 } for all i \in 1..r.
//	5) Precompute weight[i][j] = \sum_{j' >= j, f_j' \in \delta^+(e_i)} x_{f_j'}.
//	6) Get the most violated inequality with DP recursion
//		f(i, j) = max{ f(i+1, j), f(i+1, max{next(i), j}) + x_{e_i} - weight(i, j)}.
//		with boundary conditions f(i, j) = 0 if i = r + 1.
class TDFI : public goc::SeparationRoutine
{
public:
	TDFI(const VRPInstance& vrp, const goc::Matrix<std::vector<goc::Variable>>& x);
	
	virtual std::vector<goc::Constraint> Separate(const goc::Valuation& z, int node_number, int count_limit,
												  double node_bound) const;
private:
	VRPInstance vrp;
	goc::Matrix<std::vector<goc::Variable>> x;
};
} // namespace ejor2019

#endif //EJOR2019_TDFI_H
