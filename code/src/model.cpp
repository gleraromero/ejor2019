//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "model.h"

using namespace std;
using namespace goc;

namespace ejor2019
{
void Model::AddTimeWindowConstraints()
{
	for (Vertex v: vrp.D.Vertices())
	{
		f->AddConstraint((1.0*t[v]).GEQ(vrp.tw[v].left * y[v])); // Constraint (17).
		f->AddConstraint((1.0*t[v]).LEQ(vrp.tw[v].right * y[v])); // Constraint (17).
	}
}

void Model::FixStartTo0()
{
	f->SetVariableBound(t[vrp.o], 0.0, 0.0);
}

void Model::SetObjective(bool obj_duration, bool is_profitable, const vector<ProfitUnit>& profits)
{
	// Clear coefficients and set objective to minimization.
	f->Minimize(Expression());
	
	// Set new coefficients.
	if (obj_duration)
	{
		f->SetObjectiveCoefficient(t[vrp.o], -1.0);
		f->SetObjectiveCoefficient(t[vrp.d], 1.0);
	}
	if (is_profitable)
	{
		for (Vertex i: vrp.D.Vertices()) f->SetObjectiveCoefficient(y[i], -profits[i]);
	}
}

void Model::AddDurationLimitConstraint(bool fix_start_to_zero)
{
	f->AddConstraint((t[vrp.d] - t[vrp.o]).LEQ(vrp.t_max)); // Constraint (19).
	if (fix_start_to_zero)
	{
		for (Vertex i: vrp.D.Vertices())
			for (Vertex j: vrp.D.Successors(i))
				for (int m = 0; m < vrp.tau[i][j].PieceCount(); ++m)
					if (epsilon_bigger(vrp.tau[i][j].Piece(m).domain.left, vrp.t_max))
						f->SetVariableBound(xx[i][j][m], 0.0, 0.0);
	}
}

void Model::AddCapacityConstraints(bool is_pd)
{
	if (!is_pd)
	{
		f->AddConstraint(ESUM(i:vrp.D.Vertices(), vrp.q[i]*y[i]).LEQ(vrp.Q)); // Constraint (18).
	}
	else
	{
		K = vector<Variable>(vrp.D.VertexCount());
		for (Vertex i: vrp.D.Vertices()) K[i] = f->AddVariable("K_" + STR(i), VariableDomain::Real, 0.0, min(vrp.Q, vrp.Q + vrp.q[i]));
		for (Vertex i: vrp.D.Vertices())
			for (Vertex j: vrp.D.Successors(i))
				f->AddConstraint((K[i] + vrp.q[j]).LEQ(K[j] + vrp.Q * (1.0 - x[i][j]))); // Constraint (22).
	}
}

void Model::AddPickupDeliveryConstraints()
{
	int n = (vrp.D.VertexCount()-2) / 2; // Number of requests.
	for (Vertex k = 1; k <= n; ++k)
	{
		f->AddConstraint((1.0*y[k]).EQ(y[n+k])); // Constraint (20).
		f->AddConstraint((1.0*t[k]).LEQ(t[n+k])); // Constraint (21).
	}
}

void Model::ForceVisitToAllVertices()
{
	for (Vertex i: vrp.D.Vertices()) f->SetVariableBound(y[i], 1.0, 1.0);
}
} // namespace ejor2019