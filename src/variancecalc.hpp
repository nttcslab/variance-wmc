#ifndef VARIANCE_VARIANCECALC_HPP
#define VARIANCE_VARIANCECALC_HPP

#include <vector>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include "common.hpp"
#include "vtree.hpp"
#include "sdd.hpp"

/*
  VARIANCECALC class
  V: pointer to VTREE
  S: pointer to SDD
  numv: the number of vnodes in V
  mup[(literal-id)]: \mu_{P_x} value
  mun[(literal-id)]: \mu_{N_x} value
  sip[(literal-id)]: \sigma_{P_x}^2 value
  sin[(literal-id)]: \sigma_{N_x}^2 value
  exp[(node-id)]: stores E[W_{f_\alpha}]
  ev, vv, ea, va: arrays used for adjusting functions; see Algorithm 2
  cache: cache for the results of covariance computation; c[] in Algorithm 1
*/
class VARIANCECALC{
public:
  VTREE const* const V;
  SDD const* const S;
  int numv;
  std::vector<double> mup, mun, sip, sin, sipn;
  std::vector<double> exp;
  std::vector<double> ev, vv;
  std::vector<std::vector<double>> ea, va;
  std::unordered_map<std::pair<int, int>, double> cache;
  
  VARIANCECALC(VTREE const* _V, SDD const* _S): V(_V), S(_S) {
    numv = V->numv;
  };
  
  // Preprocessing; see Algorithm 2
  void preprocess();
  void preprocessInner(int now, std::vector<int>& chdrn);
  
  // ADJEXP func in Algorithm 2; ev[(vnode-id)] is stored in ea[(vnode-id)][numv]
  inline double adjustExpect(int w, int v, double val) const {
    return ea[w][v] * val;
  }
  // ADJCOV func in Algorithm 2; vv[(vnode-id)] is stored in va[(vnode-id)][numv]
  inline double adjustCovariance(int w, int v, double val, int vf, double fexp, int vg, double gexp) const {
    return va[w][v] * (val + ea[v][vf] * fexp * ea[v][vg] * gexp) + ea[w][v] * ea[w][v] * val;
  }
  
  // Compute expectation of every node; see Algorithm 3
  double computeExpect();
  double computeExpectInner(int alpha, std::vector<int>& flg);
  
  // Compute covariance; see Algorithm 1
  double compute(const std::vector<double>& _mup, const std::vector<double>& _mun, const std::vector<double>& _sip, 
                 const std::vector<double>& _sin, const std::vector<double>& _sipn);
  double computeInner(int alpha, int beta);
};

#endif  //VARIANCE_VARIANCECALC_HPP
