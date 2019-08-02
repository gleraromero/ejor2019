# A branch and cut algorithm for the time-dependent profitable tour problem with resource constraints
Source code to replicate the experiments from the article https://doi.org/10.1016/j.ejor.2019.07.014

## Abstract
In this paper we study the time-dependent profitable tour problem with resource constraints (TDPTPRC), a generalization of the profitable tour problem (PTP) which includes variable travel times to account for road congestion. In this problem, the set of customers to be served is not given and must be determined based on the profit collected when visited, keeping a balance with the total travel time. We propose a mixed integer linear programming (MILP) formulation that exploits the travel time function to reduce the size of a standard formulation from the literature. We derive four new families of valid inequalities and study the connections among them, as well as their associated separation problems. We develop a tailored Branch and Cut (BC) algorithm including these new families in addition to some well known valid inequalities from related problems. Computational results on four different problems, with alternative resources and objectives, show that the approach is flexible and effective. The algorithm achieves significant reductions in the computing times on benchmark instances from the related literature, and outperforms a recent method proposed for the time-dependent traveling salesman problem with time windows.

## Getting started
The following instructions will guide you through the steps to execute the experiments from the article.

### Prerequisites
- Python >= 3.6 [(more info)](https://www.python.org/)
- CPLEX >= 12.8 [(more info)](https://www.ibm.com/products/ilog-cplex-optimization-studio)
- Boost Graph Library >=1.66 [(more info)](https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)
    - On Linux: ```sudo apt-get install libboost-all-dev```
- CMake >= 2.8.4 [(more info)](https://cmake.org/)
    - On Linux: ```sudo apt-get install cmake```
- C++14 or higher [(more info)](https://es.wikipedia.org/wiki/C%2B%2B14)

### Built with
- Kaleidoscope: A tool to visualize the outputs of Optimization Problems [(more info)](https://github.com/gleraromero/kaleidoscope)
- Runner: A script to ease the process of running experiments [(more info)](https://github.com/gleraromero/runner)
- GOC lib: A library that includes interfaces for using (Mixed Integer) Linear Programming solvers, and some useful resources [(more info)](https://github.com/gleraromero/goc).

### Running the experiments.
1. Add environment variables with the paths to the libraries.
    1. Add two environment variables to bash with CPLEX include and library paths.
        1. ```export CPLEX_INCLUDE=<path_to_cplex_include_dir>```
            - Usually on Linux: _/opt/ibm/ILOG/CPLEX_Studio\<VERSION\>/cplex/include_
        1. ```export CPLEX_BIN=<path_to_cplex_lib_binary_file>```
            - Usually on Linux: _/opt/ibm/ILOG/CPLEX_Studio\<VERSION\>/cplex/lib/x86-64_linux/static_pic/libcplex.a_
    1. Add two environment variables to bash with BOOST Graph Library include and library paths.
        1. ```export BOOST_INCLUDE=<path_to_boost_include_dir>```
            - Usually on Linux: _/usr/include_
        1. ```export BOOST_BIN=<path_to_boost_lib_binary_file>```
            - Usually on Linux: _/usr/lib/x86_64-linux-gnu/libboost_graph.a_
1. Go to the ejor2019 root directory.
1. Execute ```python3 runner/runner.py <experiment_file>```
1. The execution output will be continually saved to the output folder.

> Experiment files are located in the _experiments_ folder. For more information see Section [Experiments](#Experiments)

### Experiments
There are eight experiment files. Which together comprise all the experiments carried out in the article.
* _Section 6.1_: tdcespp.json
* _Section 6.2_: tdop.json
* _Section 6.3_: tdptptwpd.json
* _Section 6.4_: tdtsptw.json
* _Section 6.5_: tdcespp-root.json tdop-root.json tdptptwpd-root.json tdtsptw-root.json

### Visualizing the experiment results.
1. Go to https://gleraromero.github.io/kaleidoscope/ejor2019
1. Add the output file.
1. Select the experiments.
1. Add some attributes to visualize.
1. Click on Refresh.
1. If more details on an experiment are desired click on the + icon in a specific row.

### Checker
We include a checker program to validate that algorithms produce **valid** routes. To run the checker execute:
```python3 checker/checker.py output/<output_file.json>```

The checker will go through each instance and validate:
- That the exact solution route is feasible (with respect to all resources).
- That the reported duration of the route is correct.
- If Optimum status is reported, then it should be better or equal than any solution in the _solutions.json_ file of its dataset.

## Built With
* [JSON for Modern C++](https://github.com/nlohmann/json)
* [Boost Graph Library](https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)

## Authors
- Gonzalo Lera-Romero
- Juan José Miranda-Bront

## License
This project is licensed under the MIT License - see the LICENSE.md file for details
