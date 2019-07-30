//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_capacity.h"

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace ejor2019
{
namespace
{
// Removes the arc ij from the instance.
void remove_arc(Vertex i, Vertex j, json& instance)
{
	if (instance["digraph"]["arcs"][i][j] == 0) return;
	instance["digraph"]["arcs"][i][j] = 0;
	int arc_count = instance["digraph"]["arc_count"];
	instance["digraph"]["arc_count"] = arc_count - 1;
	if (has_key(instance, "travel_times")) instance["travel_times"][i][j] = vector<json>({});
}
}

void preprocess_capacity(json& instance)
{
	double Q = instance["capacity"];
	vector<double> q = instance["demands"];
	Digraph D = instance["digraph"];
	for (Vertex i: D.Vertices())
		for (Vertex j: D.Successors(i))
			if (epsilon_bigger(q[i]+q[j], Q))
				remove_arc(i, j, instance);
}
} // namespace ejor2019