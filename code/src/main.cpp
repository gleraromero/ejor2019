//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include <iostream>
#include <vector>
#include <goc/goc.h>

#include "vrp_instance.h"
#include "preprocess_travel_times.h"
#include "preprocess_capacity.h"
#include "preprocess_time_windows.h"
#include "preprocess_service_waiting.h"

#include "model.h"
#include "ttbf.h"
#include "cttbf.h"
#include "gcs.h"
#include "weak_pi.h"
#include "weak_sigma.h"
#include "weak_pisigma.h"
#include "tdfi.h"
#include "tdcsi.h"

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace ejor2019;

// Experiment parameters.
bool is_profitable; // Indicates if the current instance is a profitable problem (i.e. vertices have profits
					// and not all have to be visited).
bool is_pd; // Indicates if the current instance is a Pickup and Delivery problem.
bool is_op; // Indicates if the current instance is an Orienteering problem.
bool is_tw; // Indicates if the current instance has time windows.
bool is_cc; // Indicates if the current instance has capacities.
bool fix_start_to_zero; // Indicates if the start time of routes should be t0==0 (for makespan minimizing instances).
bool obj_duration; 	// Indicates if minimizing duration is in the objective.
					// (In TDOP instances, objective does not include duration)
string model_type; // Model type to use "CTTBF" or "TTBF".
Duration time_limit; // Exact MIP time limit.
json solver_config; // Solver parameters.
SeparationStrategy separation_strategy; // Separation strategy chosen.

// Converts the VRP instance to a time independent one (t_ij = max(\tau_ij(t') : t' \in [0, T]) forall ij \in A).
// Solves the TI model with a small time limit and gets the BKS as the heuristic solution.
bool heuristic_solution(const VRPInstance& vrp, const vector<double>& profits,
						Route* solution, BCExecutionLog* execution_log)
{
	// Convert VRP instance to TI.
	VRPInstance vrp_ti = vrp;
	for (Arc e: vrp.D.Arcs())
	{
		Vertex i = e.tail, j = e.head;
		vrp_ti.tau[i][j] = PWLFunction();
		vrp_ti.tau[i][j].AddPiece(LinearFunction(Point2D(0.0, vrp.tau[i][j].Image().right), Point2D(vrp.T, vrp.tau[i][j].Image().right)));
	}

	// Solve model.
	CTTBF model(vrp_ti);
	model.SetObjective(obj_duration, is_profitable, profits);
	if (is_tw) model.AddTimeWindowConstraints();
	if (fix_start_to_zero) model.FixStartTo0();
	if (is_op) model.AddDurationLimitConstraint();
	if (is_cc) model.AddCapacityConstraints(is_pd);
	if (is_pd) model.AddPickupDeliveryConstraints();
	if (!is_profitable) model.ForceVisitToAllVertices();
	
	BCSolver solver;
	solver.time_limit = 5.0_sec;
	solver.config = solver_config;
	*execution_log = solver.Solve(model.f, {BCOption::RootInformation, BCOption::CutInformation, BCOption::BestIntSolution, BCOption::ScreenOutput});
	if (execution_log->best_int_solution.IsSet())
	{
		*solution = model.ParseSolution(execution_log->best_int_solution);
		solution->duration = vrp.ReadyTime(solution->path, solution->t0) - solution->t0; // Update with TD duration.
	}
	return execution_log->best_int_solution.IsSet();
}

// The experiment has the following parameters:
//	- model_type: string (formulation type)
//		* "ttbf": full formulation
//		* "cttbf": compact formulation.
//	- time_limit: number (time limit in seconds)
int main()
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_input_in_debug("instances/guerriero_et_al_2014", "15_70_A_A1", "experiments/tdtsptw.json", "CTTBF-basic (makespan)");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	is_profitable = value_or_default(experiment, "is_profitable", false);
	is_pd = value_or_default(experiment, "is_pd", false);
	is_op = value_or_default(experiment, "is_op", false);
	is_tw = value_or_default(experiment, "is_tw", false);
	is_cc = value_or_default(experiment, "is_cc", false);
	fix_start_to_zero = value_or_default(experiment, "fix_start_to_zero", false);
	obj_duration = value_or_default(experiment, "obj_duration", false);
	
	model_type = value_or_default(experiment, "model_type", "cttbf");
	time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	solver_config = value_or_default(experiment, "solver_config", json());
	separation_strategy = value_or_default(experiment, "separation_strategy", SeparationStrategy());
	
	// Show experiment details.
	clog << "Model type: " << model_type << endl;
	clog << "Time limit: " << time_limit << "sec." << endl;
	if (is_profitable) clog << "- Profitable" << endl;
	if (is_pd) clog << "- Pickup and delivery" << endl;
	if (is_op) clog << "- Orienteering" << endl;
	if (is_tw) clog << "- Has time windows" << endl;
	if (is_cc) clog << "- Has vehicle capacity" << endl;
	if (fix_start_to_zero) clog << "- Fix start to zero" << endl;
	if (obj_duration) clog << "- Minimizing duration is in the objective" << endl;
	
	// Preprocess instance JSON.
	if (is_cc) preprocess_capacity(instance);
	preprocess_travel_times(instance);
	if (is_tw) preprocess_time_windows(instance, is_profitable);
	preprocess_service_waiting(instance);
	
	// Parse instance.
	VRPInstance vrp = instance;
	
	// Read profits if present.
	vector<ProfitUnit> profits(vrp.D.VertexCount(), 0.0);
	if (is_profitable) profits = vector<ProfitUnit>(instance["profits"].begin(), instance["profits"].end());
	
	// Create model.
	Model* model = model_type == "ttbf" ? (Model*)new TTBF(vrp) : (Model*)new CTTBF(vrp);

	// Add domain specific objective and constraints.
	model->SetObjective(obj_duration, is_profitable, profits);
	if (is_tw) model->AddTimeWindowConstraints();
	if (fix_start_to_zero) model->FixStartTo0();
	if (is_op) model->AddDurationLimitConstraint();
	if (is_cc) model->AddCapacityConstraints(is_pd);
	if (is_pd) model->AddPickupDeliveryConstraints();
	if (!is_profitable) model->ForceVisitToAllVertices();

	// Create solver.
	BCSolver solver;
	solver.screen_output = &clog; // Show output in clog.
	solver.time_limit = time_limit; // Set time limit.
	solver.config = solver_config;
	solver.separation_strategy = separation_strategy;
	clog << separation_strategy << endl;
	
	// Add separation routines to the strategy.
	GCS gcs(vrp, model->x, model->y);
	WeakPi weak_pi(vrp, model->x);
	WeakSigma weak_sigma(vrp, model->x);
	WeakPiSigma weak_pisigma(vrp, model->x);
	TDFI tdfi(vrp, model->xx);
	TDCSI tdcsi(vrp, model->y, model->xx);
	solver.separation_strategy.SetSeparationRoutine("gcs", &gcs);
	solver.separation_strategy.SetSeparationRoutine("pi", &weak_pi);
	solver.separation_strategy.SetSeparationRoutine("sigma", &weak_sigma);
	solver.separation_strategy.SetSeparationRoutine("pisigma", &weak_pisigma);
	solver.separation_strategy.SetSeparationRoutine("tdfi", &tdfi);
	solver.separation_strategy.SetSeparationRoutine("tdcsi", &tdcsi);

	// Execute initial heuristic.
	clog << "Executing initial heuristic" << endl;
	Route heuristic_route;
	BCExecutionLog heuristic_log;
	if (heuristic_solution(vrp, profits, &heuristic_route, &heuristic_log))
	{
		solver.initial_solutions.push_back(model->SerializeSolution(heuristic_route));
		output["Heuristic solution"] = VRPSolution(model->f->EvaluateValuation(model->SerializeSolution(heuristic_route)), {heuristic_route});
	}
	output["Heuristic"] = heuristic_log;

	// Solve the model.
	auto log = solver.Solve(model->f, {BCOption::BestIntSolution, BCOption::CutInformation, BCOption::RootInformation});

	// Show results.
	clog << "Status: " << log.status << endl;
	if (log.best_int_solution.IsSet())
	{
		clog << "Best int: " << log.best_int_value << endl;
		Route r = model->ParseSolution(log.best_int_solution);
		clog << "Best route: " << r.path << ", departing at: " << r.t0 << ", duration: " << r.duration << endl;
		output["Best solution"] = VRPSolution(log.best_int_value, {r});
	}
	output["Exact"] = log;

	// Send JSON output to cout.
	cout << output << endl;
	return 0;
}