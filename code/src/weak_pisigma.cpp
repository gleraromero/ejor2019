//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "weak_pisigma.h"

using namespace std;
using namespace goc;

namespace ejor2019
{
WeakPiSigma::WeakPiSigma(const VRPInstance& vrp, const Matrix<Variable>& x) : vrp(vrp), x(x), P(vrp.D.VertexCount())
{
	// Initialize precedence digraph.
	P = Matrix<bool>(vrp.D.VertexCount(), vrp.D.VertexCount(), false);
	for (Vertex i: vrp.D.Vertices())
		for (Vertex j: vrp.D.Vertices())
			if (epsilon_smaller(vrp.tw[i].right, vrp.tw[j].left))
				P[i][j] = true;
	for (Vertex i: exclude(vrp.D.Vertices(), {vrp.o, vrp.d})) P[vrp.o][i] = P[i][vrp.d] = true;
}

vector<Constraint> WeakPiSigma::Separate(const Valuation& z, int node_number, int count_limit, double node_bound) const
{
	auto& D = vrp.D;
	
	// Add to N the support digraph.
	Digraph N(D.VertexCount());
	for (int i: D.Vertices())
		for (int j: D.Successors(i))
			if (epsilon_bigger(z[x[i][j]], 0.0))
				N.AddArc({i,j});
	
	// Solve max-flow with k = {1,2,...,n-1}.
	vector<Constraint> violated;
	for (Vertex i: exclude(D.Vertices(), {vrp.o, vrp.d}))
	{
		for (Vertex j: exclude(D.Vertices(), {i, vrp.o, vrp.d}))
		{
			if (!P[i][j]) continue; // We need that i < j.
			vector<bool> Q(D.VertexCount(), false); // Q[i] indicates if i \in Q.
			Q[vrp.o] = Q[vrp.d] = true;
			for (Vertex k: D.Vertices()) if (P[k][i]) Q[k] = true; // Q[k] = true if k < i.
			for (Vertex k: D.Vertices()) if (P[j][k]) Q[k] = true; // Q[k] = true if j < k.
			
			// Removing a vertex from a network equals to removing all incident arcs.
			auto c = [&] (int i, int j) { return Q[i] || Q[j] ? 0.0 : z[x[i][j]]; };
			
			double F;
			STCut ST;
			tie(F, ST) = maxflow_mincut(N, c, i, j);
			
			// Flow from i to j should be 1.0. Otherwise, there is a violated inequality.
			if (epsilon_smaller(F, 0.9))
			{
				// Calculate S - Q, V - S - Q.
				vector<Vertex> S_minus_Q, T_minus_Q;
				for (Vertex k: ST.S) if (!Q[k]) S_minus_Q.push_back(k);
				for (Vertex k: ST.T) if (!Q[k]) T_minus_Q.push_back(k);
				
				Expression left;
				for (Vertex u: S_minus_Q)
					for (Vertex v: T_minus_Q)
						if (vrp.D.IncludesArc({u,v}))
							left += 1.0 * x[u][v];
						
				violated.push_back(left.GEQ(1.0));
			}
		}
	}
	return violated;
}
} // namespace ejor2019