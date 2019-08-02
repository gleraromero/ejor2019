//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "tdfi.h"

#include <vector>

using namespace std;
using namespace goc;

namespace ejor2019
{
namespace
{
struct TDArc
{
	Vertex i, j;
	int m;
};
}

TDFI::TDFI(const VRPInstance& vrp, const Matrix<vector<Variable>>& x) : vrp(vrp), x(x)
{

}

vector<Constraint> TDFI::Separate(const Valuation& z, int node_number, int count_limit, double node_bound) const
{
	vector<Constraint> violated;
	auto DI = [&] (TDArc e) { return vrp.tau[e.i][e.j].Piece(e.m).domain; };
	auto AI = [&] (TDArc e) {
		auto& p = vrp.tau[e.i][e.j].Piece(e.m);
		return Interval(p.domain.left + p.Value(p.domain.left), p.domain.right + p.Value(p.domain.right));
	};
	
	// 	1) Fixes a vertex k \in V - {o,d}.
	for (Vertex k: exclude(vrp.D.Vertices(), {vrp.o, vrp.d}))
	{
		// 	2) Sets I = {e_1, ..., e_r}, O = {f_1, ..., f_s} the sets of arcs entering k and leaving k respectively (with z*[e_i], z*[f_i] > 0).
		vector<TDArc> I;
		for (Vertex i: vrp.D.Predecessors(k))
			for (int m = 0; m < vrp.tau[i][k].PieceCount(); ++m)
				if (epsilon_bigger(z[x[i][k][m]], 0.0))
					I.push_back(TDArc{i,k,m});
				
		vector<TDArc> O;
		for (Vertex j: vrp.D.Successors(k))
			for (int m = 0; m < vrp.tau[k][j].PieceCount(); ++m)
				if (epsilon_bigger(z[x[k][j][m]], 0.0))
					O.push_back(TDArc{k,j,m});
		
		int r = I.size(), s = O.size();
		
		//	3) Sorts I by the start of their arrival interval, and sorts O by the start of their departure interval.
		sort(I.begin(), I.end(), [&] (TDArc e1, TDArc e2) { return AI(e1).left < AI(e2).left; });
		sort(O.begin(), O.end(), [&] (TDArc e1, TDArc e2) { return DI(e1).left < DI(e2).left; });
		
		// 	4) Precompute next[i] = min { j: min(DI(f_j)) > max(AI(e_i)) or j = s + 1 } for all i \in 1..r.
		vector<int> next(r);
		for (int i = 0; i < r; ++i)
		{
			next[i] = s;
			for (int j = 0; j < s; ++j)
			{
				if (epsilon_bigger(DI(O[j]).left, AI(I[i]).right))
				{
					next[i] = j;
					break;
				}
			}
		}
		
		//	5) Precompute weight[i][j] = \sum_{j' >= j, f_j' \in \delta^+(e_i)} x_{f_j'}.
		Matrix<double> weight(r, s+1, 0.0);
		for (int i = 0; i < r; ++i)
		{
			for (int j = s-1; j >= 0; --j)
			{
				weight[i][j] = weight[i][j+1];
				if (AI(I[i]).Intersects(DI(O[j]))) // if f_j \in \delta^+(e_i).
					weight[i][j] += z[x[O[j].i][O[j].j][O[j].m]];
			}
		}
		
		//	6) Get the most violated inequality with DP recursion
		//		f(i, j) = max{ f(i+1, j), f(i+1, max{next(i), j}) + x_{e_i} - weight(i, j)}.
		//		with boundary conditions f(i, j) = 0 if i = r + 1.
		// Bottom-up DP.
		Matrix<double> f(r+1, s+1, 0.0);
		Matrix<bool> decision(r+1, s+1); // decision[i][j] = indicates if best decision is to include arc i in F.
		for (int i = r-1; i >= 0; --i)
		{
			for (int j = s - 1; j >= 0; --j)
			{
				f[i][j] = max(f[i + 1][j], f[i + 1][max(next[i], j)] + z[x[I[i].i][I[i].j][I[i].m]] - weight[i][j]);
				decision[i][j] = f[i][j] != f[i+1][j];
			}
		}
		
		// Check if there is a violated inequality.
		if (epsilon_bigger(f[0][0], 0.1))
		{
			// Reconstruct F.
			vector<TDArc> F;
			for (int i = 0, j = 0; i < r && j < s;)
			{
				if (decision[i][j]) F.push_back(I[i]);
				j = decision[i][j] ? max(next[i], j) : j;
				i = i+1;
			}
			
			// Reconstruct \delta^+(F)
			vector<TDArc> deltaF;
			for (Arc f: vrp.D.OutboundArcs(k))
			{
				for (int m = 0; m < vrp.tau[f.tail][f.head].PieceCount(); ++m)
				{
					for (TDArc e: F)
					{
						if (AI(e).Intersects(DI({f.tail, f.head, m})))
						{
							deltaF.push_back({f.tail, f.head, m});
							break;
						}
					}
				}
			}
			
			Expression left = ESUM(e: F, x[e.i][e.j][e.m]); // \sum_{e \in F} x_e
			Expression right = ESUM(f: deltaF, x[f.i][f.j][f.m]); // \sum_{f \in \delta+(F)} x_f
			
			// Validate that violation is correct.
			if (epsilon_different(left.Value(z) - right.Value(z), f[0][0]))
			{
				clog << "F" << endl;
				for (TDArc e: F) clog << e.i << " " << e.j << " " << e.m << " -> " << z[x[e.i][e.j][e.m]] << endl;
				clog << "DeltaF" << endl;
				for (TDArc e: deltaF) clog << e.i << " " << e.j << " " << e.m << " -> " << z[x[e.i][e.j][e.m]] << endl;
				
				clog << "I" << endl;
				for (TDArc e: I) clog << e.i << " " << e.j << " " << e.m << " -> " << z[x[e.i][e.j][e.m]] << endl;
				clog << "O" << endl;
				for (TDArc e: O) clog << e.i << " " << e.j << " " << e.m << " -> " << z[x[e.i][e.j][e.m]] << endl;
				clog << left << endl;
				clog << right << endl;
				clog << left.Value(z) - right.Value(z) << " vs " << f[0][0] << endl;
				
				fail("Incorrect violation TDFI");
			}
			
			violated.push_back(left.LEQ(right));
		}
	}
	return violated;
}
} // namespace ejor2019