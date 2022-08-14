#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "utility.h"
#include "objects.h"

using namespace mycode;

// utilities for generating random vector
double random_in(const double a, const double b)
{
  double r = rand() / (double)RAND_MAX;
  return (a + (b - a) * r);
}

Vector_3 random_vector_3()
{
  FT x = random_in(0.0, 1.0);
  FT y = random_in(0.0, 1.0);
  FT z = random_in(0.0, 1.0);
  return Vector_3(x, y, z);
}

Vector_2 random_vector_2()
{
  FT x = random_in(0.0, 1.0);
  FT z = random_in(0.0, 1.0);
  return Vector_2(x, z);
}
