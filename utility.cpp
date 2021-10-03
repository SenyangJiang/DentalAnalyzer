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
FT random_in(const float a, const float b)
{
  float r = rand() / (float)RAND_MAX;
  return (FT)(a + (b - a) * r);
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

// function for reading in points
void readpp(std::vector<Point_3> &points, std::string filename)
{
  std::ifstream input(filename);
  if (!input)
  {
    std::cerr << "cannot read file" << std::endl;
    return;
  }
  double x;
  double y;
  double z;
  std::string line;
  std::size_t start;
  std::size_t end;
  while (std::getline(input, line))
  {
    // std::cout << line << std::endl;
    end = line.find("<point");
    if (end != std::string::npos)
    {
      // std::cout << end << std::endl;
      while (1)
      {
        start = line.find("\"", end + 1);
        // std::cout << start << std::endl;
        if (start == std::string::npos)
        {
          break;
        }
        char var = line[start - 2];
        end = line.find("\"", start + 1);
        if (var == 'x')
        {
          x = std::stod(line.substr(start + 1, end - start - 1));
        }
        if (var == 'y')
        {
          y = std::stod(line.substr(start + 1, end - start - 1));
        }
        if (var == 'z')
        {
          z = std::stod(line.substr(start + 1, end - start - 1));
        }
      }
      // std::cout << Point_3(x,y,z) << std::endl;
      points.push_back(Point_3(x, y, z));
    }
  }
  input.close();
}
