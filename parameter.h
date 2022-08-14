#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>

using namespace std;

struct Parameter {
  string studentModel;
  string studentCenterPoint;
  string studentMidpoint;
  string studentMarginPoints;
  string studentAxialPoints;
  string studentOcclusalPoints;
  string studentGingivaPoints;
  string originalModel;
  bool divisionEnabled;
};

#endif // PARAMETER_H
