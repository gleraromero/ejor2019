//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "tdcsi.h"

#include <vector>

using namespace std;
using namespace goc;

namespace ejor2019
{

TDCSI::TDCSI(const VRPInstance& vrp, const vector<Variable>& y, const Matrix<vector<Variable>>& x) : vrp(vrp), y(y), x(x)
{
}

vector<Constraint> TDCSI::Separate(const Valuation& z, int node_number, int count_limit, double node_bound) const
{
	vector<Constraint> violated;
	
	// DI(e) = Departure interval of arc e.
	auto DI = [&] (TDArc e) { return vrp.tau[e.i][e.j].Piece(e.m).domain; };
	// AI(e) = Arrival interval of arc e.
	auto AI = [&] (TDArc e) {
		auto& p = vrp.tau[e.i][e.j].Piece(e.m);
		return Interval(p.domain.left + p.Value(p.domain.left), p.domain.right + p.Value(p.domain.right));
	};
	// Indicates if f is successor of e.
	auto is_successor = [&] (TDArc e, TDArc f) { return e.j == f.i && AI(e).Intersects(DI(f)); };
	
	//  1) Collect time-dependent arcs with z[x_e] > 0.0.
	vector<TDArc> A;
	for (Vertex i: vrp.D.Vertices())
		for (Vertex j: vrp.D.Successors(i))
			for (int m = 0; m < vrp.tau[i][j].PieceCount(); ++m)
				if (epsilon_bigger(z[x[i][j][m]], 0.0))
					A.push_back({i,j,m});
	
	// 	2) Fixes a vertex k \in V - {o,d}.
	for (Vertex k: exclude(vrp.D.Vertices(), {vrp.o, vrp.d}))
	{
		// If no flow passes through k, do not separate.
		if (epsilon_equal(z[y[k]], 0.0)) continue;
		
		// 	3) Define network D^k_flow.
		Digraph N(2+A.size()*2);
		Vertex s = N.VertexCount()-2, t = N.VertexCount()-1;
		
		// For each e \in A, we have two vertices w+(e), w-(e).
		vector<Vertex> w_plus(A.size()), w_minus(A.size());
		for (int i = 0; i < A.size(); ++i) { w_plus[i] = 2*i+1; w_minus[i] = 2*i; }
		
		// Add arcs (s, w-(ikm)) for ikm \in A.
		for (int i = 0; i < A.size(); ++i)
			if (A[i].j == k)
				N.AddArc({s, w_minus[i]});
		
		// Add arcs (w+(idm), t) for idm \in A.
		for (int i = 0; i < A.size(); ++i)
			if (A[i].j == vrp.d)
				N.AddArc({w_plus[i], t});
			
		// Add arcs (w-(e), w+(e)) for e \in A.
		for (int i = 0; i < A.size(); ++i) N.AddArc({w_minus[i], w_plus[i]});
		
		// Add arcs (w+(e), w-(f)) for e \in A, f \in \delta+(e).
		for (int i = 0; i < A.size(); ++i)
			for (int j = 0; j < A.size(); ++j)
				if (is_successor(A[i], A[j]))
					N.AddArc({w_plus[i], w_minus[j]});
		
		// c(v1, v2) = 	z(x(e)) if v1 = w+(e), v2 = w-(f) for some e \in A.
		//				(INFTY) = 2.0 otherwise.
		auto c = [&] (int v1, int v2) {
			if (v1 < s && v2 < s && v1 % 2 == 0 && v2 % 2 == 1)
			{
				TDArc e = A[v1/2];
				return z[x[e.i][e.j][e.m]];
			}
			else
			{
				return 2.0;
			}
		};
		
		double F;
		STCut ST;
		tie(F, ST) = maxflow_mincut(N, c, s, t);
		
		// flow = z(y_k) - max_violation.
		if (epsilon_bigger(z[y[k]]-F, 0.1))
		{
			// S = { A[i] : w+(i) \in C } where C a part of the C-T min-cut.
			vector<TDArc> S;
			for (int v: ST.S) if (v < s && v % 2 == 1) S.push_back(A[v/2]);
			
			// in_S[i][j][m] indicates if arc (i,j,m) \in S.
			Matrix<vector<bool>> in_S(vrp.D.VertexCount(), vrp.D.VertexCount());
			for (Arc e: vrp.D.Arcs()) in_S[e.tail][e.head] = vector<bool>(vrp.tau[e.tail][e.head].PieceCount(), false);
			for (TDArc e: S) in_S[e.i][e.j][e.m] = true;
			
			// Compute \delta+(S).
			vector<TDArc> deltaS;
			for (Vertex i: vrp.D.Vertices())
			{
				for (Vertex j: vrp.D.Successors(i))
				{
					for (int m = 0; m < vrp.tau[i][j].PieceCount(); ++m)
					{
						if (in_S[i][j][m]) continue;
						for (TDArc e: S)
						{
							if (is_successor(e, {i,j,m}))
							{
								deltaS.push_back({i,j,m});
								break;
							}
						}
					}
				}
			}
			
			// \sum_{(i,k,m) \in S \cup \delta+(S)} x_ikm <= \sum_{(i,j,m) \in \delta+(S)}
			Expression left, right;
			for (TDArc e: S) if (e.j == k) left += x[e.i][e.j][e.m];
			for (TDArc e: deltaS) if (e.j == k) left += x[e.i][e.j][e.m];
			for (TDArc e: deltaS) right += x[e.i][e.j][e.m];
			
			// Validate that violation is correct.
			if (fabs(left.Value(z) - right.Value(z) - (z[y[k]] - F)) > 0.01) fail("Incorrect violation TDCSI");
			
			violated.push_back(left.LEQ(right));
		}
	}
	return violated;
}
} // namespace ejor2019