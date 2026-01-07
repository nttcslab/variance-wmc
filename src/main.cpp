#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <chrono>
#include "common.hpp"
#include "vtree.hpp"
#include "sdd.hpp"
#include "variancecalc.hpp"

int main(int argc, char** argv){
  VTREE V;
  V.readfromFile(argv[1]);
  SDD S;
  S.readfromFile(argv[2], V);
  int numv = V.numv / 2 + 1;
  std::vector<double> mup(numv+1), mun(numv+1), sip(numv+1), sin(numv+1), sipn(numv+1);
  FILE *fp;
  if((fp = fopen(argv[3], "r")) == NULL){
    fprintf(stderr, "ERROR: reading parameter file %s failed.\n", argv[3]);
    exit(EXIT_FAILURE);
  }
  for(int i=1; i<=numv; ++i){
    fscanf(fp, "%lf%lf%lf%lf%lf", &mup[i], &mun[i], &sip[i], &sin[i], &sipn[i]);
  }
  fclose(fp);
  
  auto start = std::chrono::system_clock::now();
  V.computeLCA();
  VARIANCECALC MAN(&V, &S);
  double res = MAN.compute(mup, mun, sip, sin, sipn);
  auto finish = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count();
  printf("EXPECTATION: %.15lf\n", MAN.exp[S.root]);
  printf("VARIANCE   : %.15lf\n", res);
  fprintf(stderr, "%.lf ms\n", ctime);
  return 0;
}
