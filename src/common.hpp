#ifndef VARIANCE_COMMON_HPP
#define VARIANCE_COMMON_HPP

#include <algorithm>
#include <utility>
#include <unordered_map>
#include <cstdint>

/*
 Hash function used for the cache storing the results of covariance computation
*/
namespace std{
  template <>
  struct hash<std::pair<int, int>>{
    public:
    uint64_t operator()(const std::pair<int, int>& s) const{
      return static_cast<uint64_t>(s.first) * 314159257ULL + static_cast<uint64_t>(s.second);
    }
  };
}

#endif  //VARIANCE_COMMON_HPP
