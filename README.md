# Variance Computation for Weighted Model Counting with Knowledge Compliation Approach

This repository includes the code and data to reproduce the experimental results in the paper "Variance Computation for Weighted Model Counting with Knowledge Compliation Approach", which is accepted for AAAI 2026.
The full version of this paper will be announced upon publication on arXiv.

## How to build

All codes are written in C++.
If your environment has C++ compiler that supports (at least) C++11 features and CMake version >=3.8, you can build all binaries with the following commands:

```shell
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After that, all binaries are generated at `release/`.

## Data

All data used in the experiments are found in `data.tar.gz`.
By inflating this archive, `data/` directory appears.
The original Bayesian network descriptions are found in `data/net/`, which were retrieved from [bnRep](https://github.com/manueleleonelli/bnRep) repository.
Each network was converted into a CNF representing the encoded Boolean function by [Ace version 3.0](http://reasoning.cs.ucla.edu/ace/).
We used `-encodeOnly -noEclause -sbk05` options for Ace to ensure that the used encoding is ENC2 proposed by [Sang, Beame, and Kautz (2005)](https://aaai.org/papers/00475-aaai05-075-performing-bayesian-inference-by-weighted-model-counting/). 
All CNF files are included in `data/cnf/`, along with lmap files that describes the correspondence between the parameters and variables of Bayesian networks and the Boolean variables in the Boolean function.
Every CNF file was compiled with [SDD package](http://reasoning.cs.ucla.edu/sdd/) with no specific options other than specifying the output files.
All SDD and vtree files are included in `data/sdd/`.
The mean and variance of the weights of variables are stored in `data/params/`, which were converted from lmap files in `data/cnf/` with the `enc2lmap` program included in this repository.
The usage of `enc2lmap` is described later.

## How to reproduce experimental results

### Variance Computation Time

Variance computation can be executed by the following command:

```shell
./main [vtree_file] [sdd_file] [parameter_file]
```

`[vtree_file]`, `[sdd_file]`, and `[parameter_file]` specify the path to the file describing vtree, SDD, and parameters.
After execution, the expectation and variance of WMC of the Boolean function represented by `[sdd_file]` is computed.

_Example:_ To measure the computational time using `projectmanagement` network, run:

```shell
./main ../data/sdd/projectmanagement.net.cnf.vtree ../data/sdd/projectmanagement.net.cnf.sdd ../data/params/projectmanagement.net.cnf.param
```

Note that the variable used as a partial assignment for every network is recorded in `data/partialassignment.csv`.

### Parameter modification

The parameter files for this experiment are stored in `data/params/showcase/` directory.
For `algalactivity2` network, the parameter file for the variance computation of $\text{Pr}(\text{Chl\_a}=0)$ where the k-th parameter's variance becomes one tenth is `data/params/showcase/algalactivity2.net.cnf.param.<k>`.
Note that there are only 40 parameter files out of 43 parameters for the `algalactivity2` network because 3 parameters have value either 0 or 1 and thus they have zero variance.
For `blockchain` network, the parameter file for $\text{Pr}(\text{BA}=\text{Low})$ is `data/params/showcase/blockchain.net.cnf.param.<k>`.
For `projectmanagement` network, the parameter file for $\text{Pr}(\text{O2}=\text{YES})$ is `data/params/showcase/projectmanagement.net.cnf.O2.param.<k>` and that for $\text{Pr}(\text{O4}=\text{YES})$ is `data/params/showcase/projectmanagement.net.cnf.O4.param.<k>`.
The variance of the marginal where the k-th parameter's variance becomes one tenth can be computed by the same command.

```shell
./main [vtree_file] [sdd_file] [parameter_file]
```

_Example:_ To compute the variance of $\text{Pr}(\text{Chl\_a}=0)$ when the 7th parameter ($\text{Pr}(\text{DO}|\text{pH}=0,\text{Te}=0)$) of the `algalactivity2` network becomes one tenth, run:

```shell
./main ../data/sdd/algalactivity2.net.cnf.vtree ../data/sdd/algalactivity2.net.cnf.sdd ../data/params/showcase/algalactivity2.net.cnf.param.7
```

Note that the descriptions for parameter values is recorded in `[network_name]-parammap.csv`.

## The `enc2lmap` program

The `enc2lmap` is the program for converting an lmap file to a parameter file.
The usage is as follows:

```shell
./main [lmap_file] [theta_value] [evidence_name] <k>
```

`[lmap_file]` specifies the path to the input lmap file.
`[theta_value]` is the parameter $\theta$ of Beta distribution described in the Experiments section.
Given parameter $p$ in the data, the variance of the parameter is set to $p(1-p)/\theta$.
`[evidence_name]` specifies the name of the evidence random variable.
`<k>` is an optional argument that set the k-th parameter's variance one tenth.

## License

This software is released under the NTT license; see `LICENSE.pdf`.