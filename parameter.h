#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>

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
  bool useManualTransform;
  double transformMatrix[4][4];
};

#endif // PARAMETER_H
