#include <vector>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>
#include "common.hpp"
#include "vtree.hpp"
#include "sdd.hpp"

void SDD::readfromFile(char* filename, const VTREE& V){
  FILE *fp;
  if((fp = fopen(filename, "r")) == NULL){
    fprintf(stderr, "ERROR: reading SDD file %s failed.\n", filename);
    exit(EXIT_FAILURE);
  }
  int c, tmp;
  while((c = fgetc(fp)) != EOF){
    switch(c){
    case '\n':
      break;
    case 's':
      fscanf(fp, "dd %d", &numv);
      type.resize(numv + 1);
      rv.resize(numv + 1);
      dat.resize(numv + 1);
      ch.resize(numv + 1);
      break;
    case 'F':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 3;
      rv[tmp] = V.numv;
      fid = tmp;
      root = tmp;
      break;
    case 'T':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 2;
      rv[tmp] = V.numv;
      tid = tmp;
      root = tmp;
      break;
    case 'L':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 0;
      fscanf(fp, "%d %d", &rv[tmp], &dat[tmp]);
      root = tmp;
      break;
    case 'D':
      fscanf(fp, "%d", &tmp);
      type[tmp] = 1;
      fscanf(fp, "%d %d", &rv[tmp], &dat[tmp]);
      ch[tmp].resize(dat[tmp]);
      for(int i=0; i<dat[tmp]; ++i) fscanf(fp, "%d %d", &ch[tmp][i][0], &ch[tmp][i][1]);
      root = tmp;
      break;
    default:
      while((c = fgetc(fp)) != '\n' && c != EOF);
      break;
    }
  }
  fclose(fp);
  if(tid == -1){
    tid = numv;
    type[numv] = 2;
    rv[numv] = V.numv;
    ++numv;
  }
}

void SDD::print() const{
  for(int i=0; i<numv; ++i){
    printf("%d ", i);
    switch(type[i]){
    case 0:
      printf("L %d %d\n", rv[i], dat[i]);
      break;
    case 1:
      printf("D %d %d", rv[i], dat[i]);
      for(const auto& val : ch[i]) printf(" %d %d", val[0], val[1]);
      puts("");
      break;
    case 2:
      puts("T");
      break;
    case 3:
      puts("F");
      break;
    default:
      break;
    }
  }
  printf("root: %d\n", root);
  printf("tid : %d\n", tid);
  printf("fid : %d\n", fid);
}
