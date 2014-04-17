#pragma once

#include <cstdlib>
#include <boost/unordered_map.hpp>

struct Infection {
public:
  explicit Infection( double it, double rt ) : infT( it ), recT( rt ) {}
  double infT;
  double recT;
};

typedef boost::unordered_multimap< int, Infection > InfectionMap;
