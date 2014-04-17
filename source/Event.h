#pragma once

#include <cstdlib>
#include <set>
#include <functional>

struct Event {
 public:
     enum class Type
     {
         Death = 0,
         Fledge = 1,
         Pair = 2,
         Divorce = 3,
         Birth = 4,
         Birthday = 5,
         Infection = 6,
         Recovery = 7,
         Vaccination = 8,
         Invalid = 9,
         Last = Vaccination,
         First = Death
     };

  explicit Event(double t) : time(t), type(Type::Invalid), hostID(), s() {}
Event(double t, Type type, int hid, int stype) : time(t), type(type), hostID( hid ), s(stype) {}

  bool operator < ( const Event & rhs ) const {
    return ( time < rhs.time );
  }

  double time;
  Type type;
  int hostID;
  int s; // strain, used in infection events
};

typedef std::multiset< Event, std::less< Event > > EventPQ;
