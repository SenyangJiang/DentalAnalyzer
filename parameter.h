#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>

using namespace std;

struct Parameter {
  string studentModel;
  string studentCenterPoint;
  string studentMidpoint;
  string studentNeighborToothMarginPoint1;
  string studentNeighborToothMarginPoint2;
  string studentMarginPoints;
  string studentAxialPoints;
  string studentOcclusalPoints;
  string studentGingivaPoints;
  string originalModel;
  string originalNeighborToothMarginPoint1;
  string originalNeighborToothMarginPoint2;
  bool useManualTransform;
  bool divisionEnabled;
  double transformMatrix[4][4];
};

#endif // PARAMETER_H
