//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "weak_sigma.h"

using namespace std;
using namespace goc;

namespace ejor2019
{
WeakSigma::WeakSigma(const VRPInstance& vrp, const Matrix<Variable>& x) : vrp(vrp), x(x), P(vrp.D.VertexCount())
{
	// Initialize precedence digraph.
	P = Matrix<bool>(vrp.D.VertexCount(), vrp.D.VertexCount(), false);
	for (Vertex i: vrp.D.Vertices())
		for (Vertex j: vrp.D.Vertices())
			if (epsilon_smaller(vrp.tw[i].right, vrp.tw[j].left))
				P[i][j] = true;
	for (Vertex i: exclude(vrp.D.Vertices(), {vrp.o, vrp.d})) P[vrp.o][i] = P[i][vrp.d] = true;
}

vector<Constraint> WeakSigma::Separate(const Valuation& z, int node_number, int count_limit, double node_bound) const
{
	auto& D = vrp.D;
	
	// Add to N the support digraph.
	Digraph N(D.VertexCount());
	for (int i: D.Vertices())
		for (int j: D.Successors(i))
			if (epsilon_bigger(z[x[i][j]], 0.0))
				N.AddArc({i,j});

	// Reverse N.
	N = N.Reverse();
	
	// Solve max-flow with k = {1,2,...,n-1}.
	vector<Constraint> violated;
	for (Vertex k: exclude(D.Vertices(), {vrp.o, vrp.d}))
	{
		// Removing a vertex from a network equals to removing all incident arcs.
		auto c = [&] (int i, int j) { return P[k][i] || P[k][j] ? 0.0 : z[x[j][i]]; };
		
		double F;
		STCut ST;
		tie(F, ST) = maxflow_mincut(N, c, k, vrp.o);
		
		// Flow from i to d should be 1.0. Otherwise, there is a violated inequality.
		if (epsilon_smaller(F, 0.9))
		{
			// Calculate predecessors of S.
			vector<bool> sigma_S(vrp.D.VertexCount(), false);
			for (Vertex i: ST.S) for (Vertex j: vrp.D.Vertices()) if(P[i][j]) sigma_S[j] = true;
			
			// Calculate S - SIGMA(S), V - S - SIGMA(S).
			vector<Vertex> S_minus_sigma_S, T_minus_sigma_S;
			for (Vertex j: ST.S) if (!sigma_S[j]) S_minus_sigma_S.push_back(j);
			for (Vertex j: ST.T) if (!sigma_S[j]) T_minus_sigma_S.push_back(j);
			
			Expression left;
			for (Vertex i: T_minus_sigma_S)
				for (Vertex j: S_minus_sigma_S)
					if (vrp.D.IncludesArc({i,j}))
						left += 1.0 * x[i][j];
			
			violated.push_back(left.GEQ(1.0));
		}
	}
	return violated;
}
} // namespace ejor2019