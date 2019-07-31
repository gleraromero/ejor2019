//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_WEAK_PISIGMA_H
#define EJOR2019_WEAK_PISIGMA_H

#include <goc/goc.h>
#include "vrp_instance.h"

namespace ejor2019
{
// The pi-sigma inequalities are the following:
//	Let X, Y \subset V - {o, d} be two sets such that i < j for each i \in X, j \in Y.
//	Let Q := {o, d, \pi(X), \sigma(Y)}. Then, for any S \subset V such that X \subseteq S, Y \subseteq V - S,
//	the pi-sigma inequalities are:
//		\sum_{i \in S - Q} \sum_{j \in V - S - Q} x_ij >= 1.
//
// The weak-pi-sigma inequalities are a special case of the pi-sigma inequalities when X = {i}, Y = {j} with i < j.
// The separation routine does the following:
// 	1) Fixes two vertices i, j such that i < j.
//	2) Removes all vertices in Q = {o, d, \pi(i), \sigma(j)}.
// 	3) Solves a max-flow from i to j.
// If the flow is smaller than 1 then it means there is a violated inequality using set S from min-cut.
class WeakPiSigma : public goc::SeparationRoutine
{
public:
	WeakPiSigma(const VRPInstance& vrp, const goc::Matrix<goc::Variable>& x);
	
	virtual std::vector<goc::Constraint> Separate(const goc::Valuation& z, int node_number, int count_limit,
												  double node_bound) const;
private:
	VRPInstance vrp;
	goc::Matrix<goc::Variable> x;
	goc::Matrix<bool> P; // Precedence matrix, P[i][j] = true if i is a predecessor of j..
};
} // namespace ejor2019

#endif //EJOR2019_WEAK_PISIGMA_H
