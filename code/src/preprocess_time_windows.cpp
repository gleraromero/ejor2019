//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "preprocess_time_windows.h"

#include <vector>
#include <queue>

using namespace std;
using namespace goc;
using namespace nlohmann;

namespace ejor2019
{
namespace
{
// Calculates the time to depart to traverse arc e arriving at tf.
// Returns: INFTY if it is infeasible to depart inside the horizon.
double departing_time(const json& instance, Arc e, double tf)
{
	PWLFunction tau_e = instance["travel_times"][e.tail][e.head];
	PWLFunction arr_e = tau_e + PWLFunction::IdentityFunction(dom(tau_e));
	if (epsilon_smaller(tf, min(img(arr_e)))) return INFTY;
	else if (epsilon_bigger(tf, max(img(arr_e)))) return max(dom(arr_e));
	return arr_e.PreValue(tf);
}

// Calculates the travel time to traverse arc e departing at t0.
// Returns: INFTY if it is infeasible to arrive inside the horizon.
double travel_time(const json& instance, Arc e, double t0)
{
	PWLFunction tau_e = instance["travel_times"][e.tail][e.head];
	if (!tau_e.Domain().Includes(t0)) return INFTY;
	return tau_e(t0);
}

// Returns: the latest we can arrive to k if departing from i (and traversing arc (i, k)) without waiting.
double latest_arrival(json& instance, Vertex i, Vertex k)
{
	vector<Interval> tw = instance["time_windows"];
	if (departing_time(instance, {i, k}, tw[k].right) != INFTY) return tw[k].right;
	return tw[i].right + travel_time(instance, {i, k}, tw[i].right);
}

// Returns: the earliest we can depart from i, to reach k inside its time window without waiting.
double earliest_departure(json& instance, Vertex i, Vertex k)
{
	vector<Interval> tw = instance["time_windows"];
	if (departing_time(instance, {i, k}, tw[k].left) != INFTY)
		return departing_time(instance, {i, k}, tw[k].left) != INFTY;
	return tw[i].left;
}

// Earliest arrival time from i to all vertices if departing at a_i.
vector<double> compute_EAT(json& instance, Vertex i)
{
	return compute_earliest_arrival_time(instance["digraph"], i, instance["time_windows"][i][0], [&] (Vertex u, Vertex v, double t0) {
		return travel_time(instance, {u, v}, t0);
	});
}

// Latest departure time from all vertices to j if arriving to j at tf.
vector<double> compute_LDT(json& instance, Vertex j)
{
	return compute_latest_departure_time(instance["digraph"], j, instance["time_windows"][j][1], [&] (Vertex u, Vertex v, double t0) {
		return departing_time(instance, {u, v}, t0);
	});
}

// Removes the arc ij from the instance.
void remove_arc(json& instance, Vertex i, Vertex j)
{
	if (instance["digraph"]["arcs"][i][j] == 0) return;
	instance["digraph"]["arcs"][i][j] = 0;
	int arc_count = instance["digraph"]["arc_count"];
	instance["digraph"]["arc_count"] = arc_count - 1;
	if (has_key(instance, "travel_times")) instance["travel_times"][i][j] = vector<json>({});
}

// Returns: if the instance includes the arc.
bool includes_arc(json& instance, Arc ij)
{
	return instance["digraph"]["arcs"][ij.tail][ij.head] == 1;
}
}

void preprocess_time_windows(json& instance, bool is_profitable)
{
	Digraph D = instance["digraph"];
	int n = D.VertexCount();
	auto& V = D.Vertices();
	auto a = [&] (Vertex i) -> double { return instance["time_windows"][i][0]; };
	auto b = [&] (Vertex i) -> double { return instance["time_windows"][i][1]; };
	Vertex o = instance["start_depot"];
	Vertex d = instance["end_depot"];
	auto set_a = [&] (Vertex i, double t) { instance["time_windows"][i][0] = t; };
	auto set_b = [&] (Vertex i, double t) { instance["time_windows"][i][1] = t; };
	
	// Initialize EAT, LDT.
	Matrix<double> EAT(n,n), LDT(n,n);
	for (int i = 0; i < n; ++i) EAT[i] = compute_EAT(instance, i);
	for (int j = 0; j < n; ++j) LDT[j] = compute_LDT(instance, j);
	// Transpose LDT so LDT[i][j] is latest departure time from i to reach j.
	for (int i = 0; i < n; ++i) for (int j = i+1; j < n; ++j) swap(LDT[i][j], LDT[j][i]);
	
	// Initialize BEFORE(k) = { i | EAT(k, i) > b_i }.
	vector<vector<Vertex>> BEFORE(n);
	for (Vertex k: V)
		for (Vertex i: exclude(V, {k}))
			if (D.IncludesArc({i, k}))
				if (epsilon_bigger(EAT[k][i], b(i))) BEFORE[k].push_back(i);
	
	// Initialize AFTER(k) = { j | EAT(j, k) > b_k }.
	vector<vector<Vertex>> AFTER(n);
	for (Vertex k: V)
		for (Vertex j: exclude(V, {k}))
			if (D.IncludesArc({k, j}))
				if (epsilon_bigger(EAT[j][k], b(k))) AFTER[k].push_back(j);
	
	// Rule 1: (3.12) 	Upper bound adjustment derived from the latest arrival time at node k from its predecessors,
	//					for k \in N - {o, d}.
	for (Vertex k:exclude(V, {o,d}))
	{
		double max_arrival = -INFTY;
		for (Vertex i: D.Predecessors(k)) max_arrival = max(max_arrival, latest_arrival(instance, i, k));
		set_b(k, min(b(k), max(a(k), max_arrival)));
	}
	
	// Rule 2: (3.13)	Lower bound adjustment derived from the earliest departure time from node k to its successors,
	//					for k \in N - {o,d}.
	for (Vertex k:exclude(V, {o,d}))
	{
		double min_dep = INFTY;
		for (Vertex j: D.Successors(k)) min_dep = min(min_dep, earliest_departure(instance, k, j));
		set_a(k, max(a(k), min(b(k), min_dep)));
	}
	
	// Only for Non Profitable instances (i.e. must visit all vertices).
	if (!is_profitable)
	{
		// Rule 0: (3.11a) 	Lower bound adjustment derived from the earliest arrival time at node k from its predecessors,
		// 					for k \in N - {o,d}.
		for (Vertex k: exclude(V, {o,d}))
		{
			double max_EAT = -INFTY;
			for (Vertex i: BEFORE[k]) max_EAT = max(max_EAT, EAT[i][k]);
			set_a(k, max(a(k), max_EAT));
		}
		
		// Rule 3: (3.14a) 	Upper bound adjustment derived from the latest departure time from node k to its successors,
		// 					for k \in N - {o,d}.
		for (Vertex k: exclude(V, {o,d}))
		{
			double min_LDT = INFTY;
			for (Vertex j: AFTER[k]) min_LDT = min(min_LDT, LDT[k][j]);
			set_b(k, min(b(k), max(a(k), min_LDT)));
		}
	}
	
	// Remove infeasible tw arcs.
	for (Arc ij: D.Arcs())
	{
		int i = ij.tail, j = ij.head;
		if (epsilon_bigger(a(i)+travel_time(instance, {i, j}, a(i)), b(j))) remove_arc(instance, i, j);
	}
	
	// Only for Non Profitable instances (i.e. must visit all vertices).
	if (!is_profitable)
	{
		// Rule 4: If exists k such that EAT(k, i) > LDT(i, j) and EAT(i, j) > LDT(j, k), then remove arc (i, j).
		for (Vertex i: V)
		{
			for (Vertex j: D.Successors(i))
			{
				if (!includes_arc(instance, {i, j})) continue;
				for (Vertex k: V)
				{
					if (k == i || k == j) continue;
					if (!includes_arc(instance, {j,k})) continue;
					if (epsilon_bigger(EAT[k][i], LDT[i][j]) && epsilon_bigger(EAT[i][j], LDT[j][k]))
					{
						remove_arc(instance, i, j);
						break;
					}
				}
			}
		}
		
		// Remove transtive arcs (if i < k && k < j then remove (i, j)).
		auto predecesor = [&](Vertex i, Vertex j) {
			return i == o || j == d || epsilon_smaller(b(i), a(j)) || epsilon_bigger(EAT[j][i], b(i));
		};
		for (Vertex i: V)
			for (Vertex j: V)
				for (Vertex k: V)
					if (includes_arc(instance, {i,j}) && includes_arc(instance, {i,k}) &&
						includes_arc(instance, {k,j}) && predecesor(i, j) && predecesor(i, k) && predecesor(k, j))
						remove_arc(instance, i, j);
		
	}
}
} // namespace ejor2019