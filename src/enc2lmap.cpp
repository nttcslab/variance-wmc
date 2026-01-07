#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <utility>

int main(int argc, char **argv){
  char buf[1024];
  char nam[1024];
  
  std::vector<std::array<double, 3>> dat;
  std::vector<int> typ;
  int N;
  int enhance = -1;
  int ccnt = 0;
  int tmp1;
  double tmpd;
  
  if(argc <= 3){
    fprintf(stderr, "ERROR: too few arguments.\n");
    exit(EXIT_FAILURE);
  }
  
  double gamma = atof(argv[2]);
  std::string vname(argv[3]);
  vname += "$0";
  if(argc >= 5){
    enhance = atoi(argv[4]);
  }
  
  FILE *fp;
  if((fp = fopen(argv[1], "r")) == NULL){
    fprintf(stderr, "ERROR: loading lmap file %s failed.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  while(fgets(buf, 1023, fp) != NULL){
    if(strlen(buf) <= 3 || buf[1] != 'c') continue;
    if(buf[3] == 'N'){
      sscanf(buf, "cc$%*c$%d", &N);
      dat.resize(N+1);
      typ.resize(N+1);
    }else if(buf[3] == 'I'){
      sscanf(buf, "cc$%*c$%d$%lf$+$%s", &tmp1, &tmpd, nam);
      std::string nams(nam);
      if(nams == vname) tmpd = 0.0;
      if(tmp1 > 0){
        typ[tmp1] = 0;
        dat[tmp1][0] = tmpd;
      }else{
        dat[-tmp1][1] = tmpd;
      }
    }else if(buf[3] == 'C'){
      sscanf(buf, "cc$%*c$%d$%lf$+$", &tmp1, &tmpd);
      if(tmp1 > 0){
        typ[tmp1] = 1;
        dat[tmp1][0] = tmpd;
      }else{
        dat[-tmp1][1] = tmpd;
      }
    }
  }
  fclose(fp);
  
  for(int i=1; i<=N; ++i){
    double sig = 0.0;
    if(typ[i] == 1){
      sig = dat[i][0] * dat[i][1] / gamma;
      if(ccnt == enhance) sig /= 10.0;
      ++ccnt;
    }
    dat[i][2] = sig;
    
    printf("%.15lf %.15lf %.15lf %.15lf %.15lf\n", dat[i][0], dat[i][1], dat[i][2], dat[i][2], -dat[i][2]);
  }
  return 0;
}