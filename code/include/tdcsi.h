//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_TDCSI_H
#define EJOR2019_TDCSI_H

#include <goc/goc.h>
#include "vrp_instance.h"

namespace ejor2019
{
// The TDCSI are the following:
// Let k \in V - {o, d} be a vertex, s \subseteq {(i,j,m) \in \hat{A}, j != d}
// 		\sum_{(i,k,m) \in S \cup \delta+(S)} x_ikm <= \sum_{(i,j,m) \in \delta+(S)}
//
// The separation routine is defined in Section A2 of the Appendix.
class TDCSI : public goc::SeparationRoutine
{
public:
	TDCSI(const VRPInstance& vrp, const std::vector<goc::Variable>& y, const goc::Matrix<std::vector<goc::Variable>>& x);
	
	virtual std::vector<goc::Constraint> Separate(const goc::Valuation& z, int node_number, int count_limit,
												  double node_bound) const;
private:
	struct TDArc { int i, j, m; };
	VRPInstance vrp;
	goc::Matrix<std::vector<goc::Variable>> x;
	std::vector<goc::Variable> y;
};
} // namespace ejor2019

#endif //EJOR2019_TDCSI_H
