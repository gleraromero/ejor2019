//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_MODEL_H
#define EJOR2019_MODEL_H

#include <goc/goc.h>

#include "vrp_instance.h"

namespace ejor2019
{
// Abstraction to represent the models CTTBF and TTBF.
// Both models must implement the Model interface.
class Model
{
public:
	VRPInstance vrp;
	std::vector<goc::Variable> y; // y_i: vertex i is selected.
	goc::Matrix<goc::Variable> x; // x_ij: arc ij is selected.
	goc::Matrix<std::vector<goc::Variable>> xx; // x_ijm: arc ij is traversed during travel time zone m.
	std::vector<goc::Variable> t; // t_i: departure time from i, if y_i = 1, 0 otherwise.
	goc::Matrix<std::vector<goc::Variable>> tt; // t_ijm: t_i if x_ijm = 1, 0 otherwise.
	std::vector<goc::Variable> K; // K_i: vehicle load at vertex i (for PD variant) (variable Q_i in the article).
	goc::Formulation* f;
	
	// Adds time window constraints (17).
	void AddTimeWindowConstraints();
	
	// For makespan objective: fixes t_o = 0.
	// Explained after constraint (19) in Section 3.3.
	void FixStartTo0();
	
	// Sets the objective function.
	// 	- obj_duration: indicates if minimizing the duration is included in the objective.
	//	- is_profitable: indicates if maximizing the profit is included in the objective.
	//	- profits: in case is_profitable == true then the profits for the vertices.
	void SetObjective(bool obj_duration, bool is_profitable, const std::vector<ProfitUnit>& profits = {});
	
	// Adds the duration limit constraint (19).
	// If the start is fixed to 0, then set all variables t_ijm = 0 with w_ijm > t_max.
	void AddDurationLimitConstraint(bool fix_start_to_zero);
	
	// Adds the capacity constraints (18) if is_pd = 0, (22)-(23) if is_pd = 1.
	//	- is_pd: indicates if the problem is a Pickup and Delivery problem.
	void AddCapacityConstraints(bool is_pd);
	
	// Adds the one-to-one pickup and delivery constraints (20)-(21).
	void AddPickupDeliveryConstraints();
	
	// If the instance is a TSP then all vertices must be visited exactly once.
	void ForceVisitToAllVertices();
	
	// Serializes the route into a valuation for the formulation.
	virtual goc::Valuation SerializeSolution(const goc::Route& route) const = 0;
	
	// Parses a solution from a solver.
	// Returns: the solution route.
	virtual goc::Route ParseSolution(const goc::Valuation& z) const = 0;
};
} // namespace ejor2019

#endif //EJOR2019_MODEL_H
