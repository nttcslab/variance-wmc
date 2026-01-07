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

void VTREE::readfromFile(char* filename){
  FILE *fp;
  if((fp = fopen(filename, "r")) == NULL){
    fprintf(stderr, "ERROR: reading vtree file %s failed.\n", filename);
    exit(EXIT_FAILURE);
  }
  int c, tmp;
  while((c = fgetc(fp)) != EOF){
    switch(c){
    case '\n':
      break;
    case 'v':
      fscanf(fp, "tree %d", &numv);
      type.resize(numv);
      dat.resize(numv);
      par.resize(numv);
      std::fill(par.begin(), par.end(), -1);
      break;
    case 'L':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 0;
      fscanf(fp, "%d", &dat[tmp][0]);
      root = tmp;
      break;
    case 'I':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 1;
      fscanf(fp, "%d %d", &dat[tmp][0], &dat[tmp][1]);
      par[dat[tmp][0]] = tmp;
      par[dat[tmp][1]] = tmp;
      root = tmp;
    default:
      while((c = fgetc(fp)) != '\n' && c != EOF);
      break;
    }
  }
  fclose(fp);
}

void VTREE::print() const{
  for(int i=0; i<numv; ++i){
    printf("%d %c %d ", i, type[i] ? 'I' : 'L', par[i]);
    if(type[i] == 0) printf("%d\n", dat[i][0]);
    else             printf("%d %d\n", dat[i][0], dat[i][1]);
  }
  printf("root: %d\n", root);
}

/*
Instead of using sophisticated LCA algorithms that uses linear preprocessing time, we store the LCA of any pair of vnodes explicitly;
although that costs quadratic time in the size of vtree, it is absorbed in the preprocessing time.
*/
void VTREE::computeLCA(){
  lca.resize(numv+1);
  for(auto&& vec : lca) vec.resize(numv+1);
  std::vector<int> chdrn;
  chdrn.reserve(numv);
  computeLCAInner(root, chdrn);
  for(int i=0; i<numv; ++i) lca[i][numv] = lca[numv][i] = i;
  lca[numv][numv] = numv;
}

void VTREE::computeLCAInner(int now, std::vector<int>& chdrn){
  lca[now][now] = now;
  if(type[now] == 0){
    chdrn.emplace_back(now);
    return;
  }
  std::vector<int> chdrnl, chdrnr;
  chdrnl.reserve(numv);
  chdrnr.reserve(numv);
  computeLCAInner(dat[now][0], chdrnl);
  computeLCAInner(dat[now][1], chdrnr);
  for(const auto& val : chdrnl)for(const auto& val2 : chdrnr) lca[val][val2] = lca[val2][val] = now;
  for(const auto& val : chdrnl){
    lca[val][now] = lca[now][val] = now;
    chdrn.emplace_back(val);
  }
  for(const auto& val : chdrnr){
    lca[val][now] = lca[now][val] = now;
    chdrn.emplace_back(val);
  }
  chdrn.emplace_back(now);
  return;
}
