//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_WEAK_PI_H
#define EJOR2019_WEAK_PI_H

#include <goc/goc.h>
#include "vrp_instance.h"

namespace ejor2019
{
// The weak-pi inequalities are the following:
//	Given a vertex k \subset V - {o,d}, k \in S \subset V - {o,d}.
// 	\sum_{i \in S - PI(S)} \sum_{j \in V - S - PI(S)} x_ij >= 1.
//
// The separation routine does the following:
// 	1) Fixes a vertex i.
// 	2) Removes all vertices k \in PI(i).
// 	3) Solves a max-flow from i to d.
// If the flow is smaller than 1 then it means there is a violated inequality using set S from min-cut.
class WeakPi : public goc::SeparationRoutine
{
public:
	WeakPi(const VRPInstance& vrp, const goc::Matrix<goc::Variable>& x);
	
	virtual std::vector<goc::Constraint> Separate(const goc::Valuation& z, int node_number, int count_limit,
		double node_bound) const;
private:
	VRPInstance vrp;
	goc::Matrix<goc::Variable> x;
	goc::Matrix<bool> P; // Precedence matrix, P[i][j] = true if i is a predecessor of j..
};
} // namespace ejor2019

#endif //EJOR2019_WEAK_PI_H
