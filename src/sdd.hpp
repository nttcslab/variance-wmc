#ifndef VARIANCE_SDD_HPP
#define VARIANCE_SDD_HPP

#include <vector>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>
#include "common.hpp"
#include "vtree.hpp"

/*
  SDD class: Node id starts at 0, while literal id starts at 1
  numv: the number of nodes in SDD
  root: the id of root node
  tid : the id of true node
  fid : the id of false node
  type[(node-id)]: 0 for literal, 1 for internal (decomp.) node, 2 for true, 3 for false
  rv[(node-id)]  : decomposition vnode id
  dat[(node-id)] : for literal node, store the literal id; multiplied by (-1) if it is a negative literal
                   for internal node, store the number of prime-sub pairs
  ch[(node-id)]  : for internal node, store the ids of child nodes of prime-sub pairs
  
  Assuming that the decomposition vnode id for true and false nodes is #(vnodes)
*/
class SDD{
public:
  int numv, root;
  int tid = -1, fid = -1;
  std::vector<int> type, rv, dat;
  std::vector<std::vector<std::array<int, 2>>> ch;
  
  // read SDD file output from SDD package
  void readfromFile(char* filename, const VTREE& V);
  // print SDD
  void print() const;
};

#endif  //VARIANCE_SDD_HPP
