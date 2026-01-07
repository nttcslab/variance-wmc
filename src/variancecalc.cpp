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
#include "variancecalc.hpp"

/*
While Algorithm 2 in the paper proceeds in a bottom-up manner along the vtree, the implementation here proceeds
in a top-down manner along the vtree; both approaches correctly work.
*/
void VARIANCECALC::preprocess(){
  ev.resize(numv);
  vv.resize(numv);
  ea.resize(numv+1);
  for(auto&& vec : ea) vec.resize(numv+1);
  va.resize(numv+1);
  for(auto&& vec : va) vec.resize(numv+1);
  std::vector<int> chdrn;
  chdrn.reserve(numv);
  preprocessInner(V->root, chdrn);
  for(int i=0; i<numv; ++i){
    ea[i][numv] = ev[i];
    va[i][numv] = vv[i];
  }
  ea[numv][numv] = 1.0;
  va[numv][numv] = 0.0;
}

void VARIANCECALC::preprocessInner(int now, std::vector<int>& chdrn){
  ea[now][now] = 1;
  va[now][now] = 0;
  if(V->type[now] == 0){
    int vleaf = V->dat[now][0];
    chdrn.emplace_back(now);
    ev[now] = mup[vleaf] + mun[vleaf];
    vv[now] = sip[vleaf] + sin[vleaf] + sipn[vleaf] * 2.0;
    return;
  }
  int idl = V->dat[now][0], idr = V->dat[now][1];
  std::vector<int> chdrnl, chdrnr;
  chdrnl.reserve(numv);
  chdrnr.reserve(numv);
  preprocessInner(idl, chdrnl);
  preprocessInner(idr, chdrnr);
  ev[now] = ev[idl] * ev[idr];
  vv[now] = vv[idl] * vv[idr] + vv[idl] * ev[idr] * ev[idr] + ev[idl] * ev[idl] * vv[idr];
  for(const auto& val : chdrnl){
    ea[now][val] = ea[idl][val] * ev[idr];
    va[now][val] = va[idl][val] * (vv[idr] + ev[idr] * ev[idr]) + ea[idl][val] * ea[idl][val] * vv[idr];
    chdrn.emplace_back(val);
  }
  for(const auto& val : chdrnr){
    ea[now][val] = ev[idl] * ea[idr][val];
    va[now][val] = vv[idl] * (va[idr][val] + ea[idr][val] * ea[idr][val]) + ev[idl] * ev[idl] * va[idr][val];
    chdrn.emplace_back(val);
  }
  chdrn.emplace_back(now);
  return;
}

/*
See Algorithm 3; exp[] corresponds to e[] in Algorithm 3
*/
double VARIANCECALC::computeExpect(){
  std::vector<int> flg(S->numv);
  exp.reserve(S->numv);
  exp[S->tid] = 1.0;
  flg[S->tid] = 1;
  if(S->fid != -1){
    exp[S->fid] = 0.0;
    flg[S->fid] = 1;
  }
  return computeExpectInner(S->root, flg);
}

double VARIANCECALC::computeExpectInner(int alpha, std::vector<int>& flg){
  // cache hit (line 1)
  if(flg[alpha]) return exp[alpha];
  // base cases (lines 3--5)
  if(S->type[alpha] == 2) exp[alpha] = 1;
  else if(S->type[alpha] == 3) exp[alpha] = 0;
  else if(S->type[alpha] == 0) exp[alpha] = S->dat[alpha] > 0 ? mup[S->dat[alpha]] : mun[-S->dat[alpha]];
  // In SDD, \vee and \wedge-nodes must alternate. We use this property to simplify the implementation.
  else{
    exp[alpha] = 0;
    for(const auto& pr : S->ch[alpha]){
      double tmpl = computeExpectInner(pr[0], flg);
      double tmpr = computeExpectInner(pr[1], flg);
      exp[alpha] += adjustExpect(V->dat[S->rv[alpha]][0], S->rv[pr[0]], tmpl) * adjustExpect(V->dat[S->rv[alpha]][1], S->rv[pr[1]], tmpr);
    }
  }
  flg[alpha] = 1;
  return exp[alpha];
}

/*
See Algorithm 1; V[W_f] can be computed by Cov[W_f,W_f]
*/
double VARIANCECALC::compute(const std::vector<double>& _mup, const std::vector<double>& _mun, const std::vector<double>& _sip, 
                             const std::vector<double>& _sin, const std::vector<double>& _sipn){
  mup = _mup;
  mun = _mun;
  sip = _sip;
  sin = _sin;
  sipn = _sipn;
  preprocess();
  computeExpect();
  cache.reserve(S->numv * S->numv);
  return computeInner(S->root, S->root);
}

double VARIANCECALC::computeInner(int alpha, int beta){
  // base cases (lines 3--4)
  if(alpha == S->fid || beta == S->fid) return 0;
  if(alpha == S->tid && beta == S->tid) return 0;
  // search for cache (line 1)
  std::pair<int, int> key(std::min(alpha, beta), std::max(alpha, beta));
  auto itr = cache.find(key);
  if(itr != cache.end()){
    return itr->second;
  }
  int da = S->rv[alpha];
  int db = S->rv[beta];
  int anc = V->lca[da][db];
  int typa, typb;
  typa = (da == numv) ? 3 : ((anc == da) ? 0 : ((V->lca[V->dat[anc][0]][da] == anc) ? 2 : 1));
  typb = (db == numv) ? 3 : ((anc == db) ? 0 : ((V->lca[V->dat[anc][0]][db] == anc) ? 2 : 1));
  double res = 0.0;
  // literal cases (lines 11--16)
  if(V->type[anc] == 0){
    int tmpa, tmpb;
    if(da == numv){
      tmpb = S->dat[beta];
      res = (tmpb > 0) ? (sip[tmpb] + sipn[tmpb]) : (sin[-tmpb] + sipn[-tmpb]);
    }else if(db == numv){
      tmpa = S->dat[alpha];
      res = (tmpa > 0) ? (sip[tmpa] + sipn[tmpa]) : (sin[-tmpa] + sipn[-tmpa]);
    }else{
      tmpa = S->dat[alpha];
      tmpb = S->dat[beta];
      if     (tmpa > 0 && tmpb > 0) res = sip[tmpa];
      else if(tmpa < 0 && tmpb < 0) res = sin[-tmpa];
      else                          res = sipn[std::max(tmpa, tmpb)];
    }
  }
  // We simplify the implementation by using the property of SDD that \vee and \wedge-nodes must alternate
  else{
    int ancl = V->dat[anc][0];
    int ancr = V->dat[anc][1];
    std::vector<std::array<int, 2>> const* cha;
    std::vector<std::array<int, 2>> const* chb;
    std::vector<std::array<int, 2>> tmpva(1);
    std::vector<std::array<int, 2>> tmpvb(1);
    switch(typa){
    case 0: cha = &(S->ch[alpha]); break;
    case 1: tmpva[0][0] = alpha; tmpva[0][1] = S->tid; cha = &tmpva; break;
    case 2: tmpva[0][0] = S->tid; tmpva[0][1] = alpha; cha = &tmpva; break;
    default: tmpva[0][0] = tmpva[0][1] = S->tid; cha = &tmpva; break;
    }
    switch(typb){
    case 0: chb = &(S->ch[beta]); break;
    case 1: tmpvb[0][0] = beta; tmpvb[0][1] = S->tid; chb = &tmpvb; break;
    case 2: tmpvb[0][0] = S->tid; tmpvb[0][1] = beta; chb = &tmpvb; break;
    default: tmpvb[0][0] = tmpvb[0][1] = S->tid; chb = &tmpvb; break;
    }
    double restmp = 0;
    for(const auto& pa : *cha)for(const auto& pb : *chb){
      if(pa[1] == S->fid || pb[1] == S->fid) restmp = 0;
      else{
        int rval = S->rv[pa[0]], rvbl = S->rv[pb[0]], rvar = S->rv[pa[1]], rvbr = S->rv[pb[1]];
        double el = adjustExpect(ancl, rval, exp[pa[0]]) * adjustExpect(ancl, rvbl, exp[pb[0]]);
        double er = adjustExpect(ancr, rvar, exp[pa[1]]) * adjustExpect(ancr, rvbr, exp[pb[1]]);
        double resl = computeInner(pa[0], pb[0]);
        double resr = computeInner(pa[1], pb[1]);
        double cl = adjustCovariance(ancl, V->lca[rval][rvbl], resl, rval, exp[pa[0]], rvbl, exp[pb[0]]);
        double cr = adjustCovariance(ancr, V->lca[rvar][rvbr], resr, rvar, exp[pa[1]], rvbr, exp[pb[1]]);
        restmp = cl * (cr + er) + el * cr;
      }
      res += restmp;
    }
  }
  cache.emplace(key, res);
  return res;
}