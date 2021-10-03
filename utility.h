#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include "objects.h"

// utilities for generating random vector
mycode::FT random_in(const float a, const float b);

mycode::Vector_3 random_vector_3();

mycode::Vector_2 random_vector_2();

// read in points in .pp file
void readpp(std::vector<mycode::Point_3> &points, std::string filename);

#endif
