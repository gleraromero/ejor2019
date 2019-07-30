//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_CTTBF_H
#define EJOR2019_CTTBF_H

#include <vector>

#include <goc/goc.h>

#include "model.h"
#include "vrp_instance.h"

namespace ejor2019
{
// Model TTBF-Compact from the article.
class CTTBF : public Model
{
public:
	goc::Matrix<goc::Variable> t_hat; // t_hat_ij: t_i, if x_ijm = 1 with theta(ijm)=0.
	
	// Initializes a formulation f with the specified parameters.
	CTTBF(const VRPInstance& vrp);
	
	~CTTBF();
	
	// Serializes the route into a valuation for the formulation.
	virtual goc::Valuation SerializeSolution(const goc::Route& route) const;
	
	// Parses a solution from a solver.
	// Returns: the solution route.
	virtual goc::Route ParseSolution(const goc::Valuation& z) const;
};
} // namespace ejor2019
#endif //EJOR2019_CTTBF_H
