//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "cttbf.h"

using namespace std;
using namespace goc;

namespace ejor2019
{
CTTBF::CTTBF(const VRPInstance& vrp)
{
	this->vrp = vrp;
	int n = vrp.D.VertexCount();
	auto& V = vrp.D.Vertices(); auto& A = vrp.D.Arcs();
	auto in = [&] (Vertex j) { return vrp.D.Predecessors(j); };
	auto out = [&] (Vertex i) { return vrp.D.Successors(i); };
	auto theta = [&] (Vertex i, Vertex j, int m) { return vrp.tau[i][j].Piece(m).slope; };
	auto eta = [&] (Vertex i, Vertex j, int m) { return vrp.tau[i][j].Piece(m).intercept; };
	auto w = [&] (Vertex i, Vertex j, int m) { return m == 0 ? vrp.tau[i][j].Piece(0).domain.left : vrp.tau[i][j].Piece(m-1).domain.right; };
	auto& q = vrp.q; auto Q = vrp.Q;
	Vertex o = vrp.o, d = vrp.d; // origin and destination.
	vector<Vertex> Vo = exclude(V, {o}); // Vo = V - {o}.
	vector<Vertex> Vd = exclude(V, {d}); // Vd = V - {d}.
	vector<Vertex> Vod = exclude(V, {o, d}); // Vod = V - {o, d}.
	Matrix<vector<int>> M(n, n); // M[i][j] = { m \in 0...|T^ij_m|-1 }
	for (Arc e: A) M[e.tail][e.head] = range(0, vrp.tau[e.tail][e.head].PieceCount());
	Matrix<vector<int>> M_nz(n, n); // M_nz[i][j] = { m \in M : theta(ijm) != 0 }
	for (auto i: V) for (auto j: out(i)) for (int m: M[i][j]) if (epsilon_different(theta(i,j,m),0.0)) M_nz[i][j].push_back(m);
	Matrix<vector<int>> M_z(n, n); // M_nz[i][j] = { m \in M : theta(ijm) == 0 }
	for (auto i: V) for (auto j: out(i)) for (int m: M[i][j]) if (epsilon_equal(theta(i,j,m),0.0)) M_z[i][j].push_back(m);
	int M_max = 0;
	for (Arc e: A) M_max = max(M_max, (int)M[e.tail][e.head].size());
	
	// Create formulation.
	f = BCSolver::NewFormulation();
	
	// Add variables to formulation.
	y = vector<Variable>(n);
	x = Matrix<Variable>(n, n);
	xx = Matrix<vector<Variable>>(n, n, vector<Variable>(M_max));
	t = vector<Variable>(n);
	t_hat = Matrix<Variable>(n,n);
	tt = Matrix<vector<Variable>>(n, n, vector<Variable>(M_max));
	for (Vertex i: V) y[i] = f->AddVariable("y_" + STR(i), VariableDomain::Binary, 0.0, 1.0);
	for (Arc e: A) x[e.tail][e.head] = f->AddVariable("x_" + STR(e.tail) + "_" + STR(e.head), VariableDomain::Binary, 0.0, 1.0);
	for (Arc e: A) for (int m: M[e.tail][e.head]) xx[e.tail][e.head][m] = f->AddVariable("x_" + STR(e.tail) + "_" + STR(e.head) + "_" + STR(m), VariableDomain::Binary, 0.0, 1.0);
	for (Vertex i: V) t[i] = f->AddVariable("t_" + STR(i), VariableDomain::Real, 0.0, vrp.T);
	for (auto i: V) for (auto j: out(i)) t_hat[i][j] = f->AddVariable("th_" + STR(i) + "_" + STR(j), VariableDomain::Real, 0.0, vrp.T);
	for (Arc e: A) for (int m: M_nz[e.tail][e.head]) tt[e.tail][e.head][m] = f->AddVariable("t_" + STR(e.tail) + "_" + STR(e.head) + "_" + STR(m), VariableDomain::Real, 0.0, vrp.T);
	
	// Constraints.
	for (auto i: V) for (auto j: out(i)) f->AddConstraint(ESUM(m: M[i][j], xx[i][j][m]).EQ(x[i][j])); // (2)
	for (auto j: Vo) f->AddConstraint(ESUM(i:in(j), x[i][j]).EQ(y[j])); // (3)
	f->AddConstraint(ESUM(j:out(o), x[o][j]).EQ(ESUM(i:in(d), x[i][d]))); // (4)
	f->AddConstraint(ESUM(j:out(o), x[o][j]).EQ(1.0)); // (4)
	for (auto k: Vod) f->AddConstraint((ESUM(i:in(k), x[i][k]) - ESUM(j:out(k), x[k][j])).EQ(0)); // (5)
	for (auto j: Vo) f->AddConstraint(ESUM(i:in(j), t_hat[i][j] + ESUM(m:M_nz[i][j], (1+theta(i,j,m)) * tt[i][j][m]) + ESUM(m:M[i][j], eta(i,j,m) * xx[i][j][m])).LEQ(t[j])); // (6)
	for (auto i: Vd) f->AddConstraint(ESUM(j:out(i), t_hat[i][j] + ESUM(m:M_nz[i][j], tt[i][j][m])).EQ(t[i])); // (7)
	for (auto i: V) for (auto j: out(i)) for(int m: M_nz[i][j]) f->AddConstraint((w(i,j,m) * xx[i][j][m]).LEQ(tt[i][j][m])); // (8)
	for (auto i: V) for (auto j: out(i)) for(int m: M_nz[i][j]) f->AddConstraint((w(i,j,m+1) * xx[i][j][m]).GEQ(tt[i][j][m])); // (8)
	for (auto i: V) for (auto j: out(i)) f->AddConstraint(ESUM(m:M_z[i][j], w(i,j,m) * xx[i][j][m]).LEQ(t_hat[i][j])); // (8)
	for (auto i: V) for (auto j: out(i)) f->AddConstraint(ESUM(m:M_z[i][j], w(i,j,m+1) * xx[i][j][m]).GEQ(t_hat[i][j])); // (8)
	
	// GCS with |S|=2.
	for (auto i: V) for (auto j: out(i)) if (vrp.D.IncludesArc({j,i})) f->AddConstraint((x[i][j]+x[j][i]).LEQ(y[i]));
	for (auto i: V) for (auto j: out(i)) if (vrp.D.IncludesArc({j,i})) f->AddConstraint((x[i][j]+x[j][i]).LEQ(y[j]));
}

CTTBF::~CTTBF()
{
	delete f;
}

Valuation CTTBF::SerializeSolution(const Route& r) const
{
	const GraphPath& p = r.path;
	TimeUnit t0 = r.t0;
	
	Valuation z;
	for (auto v: p) z.SetValue(y[v], 1.0);
	for (int k = 0; k < (int)p.size()-1; ++k) z.SetValue(x[p[k]][p[k+1]], 1.0);
	for (int k = 0; k < (int)p.size()-1; ++k)
	{
		Vertex i = p[k], j = p[k+1];
		z.SetValue(t[i], t0);
		int m = vrp.tau[i][j].PieceIncluding(t0);
		if (epsilon_equal(vrp.tau[i][j].Piece(m).slope, 0.0)) z.SetValue(t_hat[i][j], t0);
		else z.SetValue(tt[i][j][m], t0);
		z.SetValue(xx[i][j][m], 1.0);
		t0 += vrp.tau[i][j](t0);
	}
	z.SetValue(t[p.back()], t0);
	return z;
}

Route CTTBF::ParseSolution(const Valuation& z) const
{
	double t0 = z[t[vrp.o]];
	GraphPath p = {vrp.o};
	while (p.back() != vrp.d)
	{
		Vertex u = p.back();
		for (auto v: vrp.D.Successors(u))
			if (epsilon_equal(z[x[u][v]], 1.0))
				p.push_back(v);
	}
	return Route(p, t0, vrp.ReadyTime(p, t0)-t0);
}
} // namespace euroalio