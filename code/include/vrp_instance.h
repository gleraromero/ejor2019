//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef EJOR2019_VRP_INSTANCE_H
#define EJOR2019_VRP_INSTANCE_H

#include <vector>
#include <goc/goc.h>

namespace ejor2019
{
typedef double TimeUnit; // Represents time.
typedef double CapacityUnit; // Represents the capacity.
typedef double ProfitUnit; // Represents the profit of vertices.

// This class represents an instance of a vehicle routing problem.
// Considerations:
// 	- It considers two depots (origin and destination).
class VRPInstance : public goc::Printable
{
public:
	goc::Digraph D; // digraph representing the network.
	goc::Vertex o, d; // origin and destination depot.
	TimeUnit T; // end of planning horizon ([0,T]).
	std::vector<goc::Interval> tw; // time window of customers (tw[i] = time window of customer i).
	CapacityUnit Q; // vehicle capacity.
	std::vector<CapacityUnit> q; // demand of customers (q[i] = demand of customer i).
	TimeUnit t_max; // maximum duration of a path.
	goc::Matrix<goc::PWLFunction> tau; // tau[i][j](t) = travel time of arc (i, j) if departing from i at t.
	
	// Returns: the time we finish visiting the last vertex if departing at t0.
	// If infeasible, returns INFTY.
	TimeUnit ReadyTime(const goc::GraphPath& p, TimeUnit t0=0) const;
	
	// Prints the JSON representation of the instance.
	virtual void Print(std::ostream& os) const;
};

// Serializes the instance.
void to_json(nlohmann::json& j, const VRPInstance& instance);

// Parses an instance.
void from_json(const nlohmann::json& j, VRPInstance& instance);
} // namespace ejor2019

#endif //EJOR2019_VRP_INSTANCE_H
