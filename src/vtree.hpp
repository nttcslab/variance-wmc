#ifndef VARIANCE_VTREE_HPP
#define VARIANCE_VTREE_HPP

#include <vector>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include "common.hpp"

/*
  VTREE class: Vnode id starts at 0, while literal id starts at 1
  numv: the number of vnodes in vtree
  root: the id of root vnode
  type[(vnode-id)]: 0 for leaf vnode, 1 for internal vnode
  dat[(vnode-id)]: for leaf vnode, dat[][0] stores literal id;
                   for internal vnode, dat[][0] and dat[][1] stores ids of child vnodes
  lca[(vnode-id1)][(vnode-id2)]: the id of the lowest common ancestor of two vnodes
*/
class VTREE{
public:
  int numv, root;
  std::vector<int> type;
  std::vector<int> par;
  std::vector<std::array<int, 2>> dat;
  std::vector<std::vector<int>> lca;
  
  // read vtree file output from SDD package
  void readfromFile(char* filename);
  // print vtree
  void print() const;
  // compute LCA for any pair of vnodes as a preprocessing
  void computeLCA();
  void computeLCAInner(int now, std::vector<int>& chdrn);
};

#endif  //VARIANCE_VTREE_HPP
