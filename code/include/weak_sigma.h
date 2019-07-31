//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_WEAK_SIGMA_H
#define EJOR2019_WEAK_SIGMA_H

#include <goc/goc.h>
#include "vrp_instance.h"

namespace ejor2019
{
// The sigma inequalities are the following:
//	Let S \subset V - {o, d} be a set of vertices, then the sigma-inequalities are:
// 		\sum_{i \in V - S - \sigma(S)} \sum_{j \in S - \sigma(S)} x_ij >= 1.
//
// The weak-sigma inequalities are a special case of the sigma inequalities.
// The separation routine does the following:
// 	1) Fixes a vertex i.
//	2) Reverses the support digraph.
// 	3) Removes all vertices k \in \sigma(i).
// 	4) Solves a max-flow from i to o.
// If the flow is smaller than 1 then it means there is a violated inequality using set S from min-cut.
class WeakSigma : public goc::SeparationRoutine
{
public:
	WeakSigma(const VRPInstance& vrp, const goc::Matrix<goc::Variable>& x);
	
	virtual std::vector<goc::Constraint> Separate(const goc::Valuation& z, int node_number, int count_limit,
												  double node_bound) const;
private:
	VRPInstance vrp;
	goc::Matrix<goc::Variable> x;
	goc::Matrix<bool> P; // Precedence matrix, P[i][j] = true if i is a predecessor of j..
};
} // namespace ejor2019

#endif //EJOR2019_WEAK_SIGMA_H
