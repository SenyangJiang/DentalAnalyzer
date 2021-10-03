#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/pointmatcher/register_point_sets.h>
#include <CGAL/pointmatcher/compute_registration_transformation.h>
#include <CGAL/OpenGR/compute_registration_transformation.h>

#include <CGAL/Polygon_mesh_processing/transform.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <unordered_set>
#include <stdlib.h>
#include <time.h>

#include <CGAL/intersections.h>

#include "utility.h"
#include "compare.h"
#include "objects.h"

using namespace mycode;

// main function
int main(int argc, const char **argv)
{
  srand(time(NULL)); // generate random seed
  std::vector<Pwn> pwns1, pwns2;

  const char *fname1 = "data/AndreCataluna.off";
  const char *fname2 = "data/KimChoi.off";
  const char *margin_file_1 = "data/AndreCataluna_margin_points.pp";
  const char *margin_file_2 = "data/KimChoi_margin_points.pp";

  std::unordered_set<vertex_descriptor> vertexSet1;
  std::unordered_set<vertex_descriptor> vertexSet2;

  select_tooth_points(fname1, margin_file_1, pwns1, vertexSet1);
  select_tooth_points(fname2, margin_file_2, pwns2, vertexSet2);

  std::cerr << "Computing registration transformation using OpenGR Super4PCS.." << std::endl;
  // First, compute registration transformation using OpenGR Super4PCS
  K::Aff_transformation_3 res1 =
      std::get<0>( // get first of pair, which is the transformation
          CGAL::OpenGR::compute_registration_transformation(pwns1, pwns2,
                                                            params::point_map(Point_map()).normal_map(Normal_map()),
                                                            params::point_map(Point_map()).normal_map(Normal_map())));
  std::cerr << "Computing registration transformation using PointMatcher ICP, "
            << "taking transformation computed by OpenGR Super4PCS as initial transformation.." << std::endl;
  K::Aff_transformation_3 res2 =
      std::get<0>(
          CGAL::pointmatcher::compute_registration_transformation(pwns1, pwns2,
                                                                  params::point_map(Point_map()).normal_map(Normal_map()),
                                                                  params::point_map(Point_map()).normal_map(Normal_map()).transformation(res1) /* initial transform for pwns2.
  				       * a transform returned from a coarse registration algorithm.
  				       * */
                                                                  ));

  // read in master model
  std::ifstream input(fname1);
  if (!input)
  {
    std::cerr << "unable to open file" << std::endl;
    return EXIT_FAILURE;
  }
  Polyhedron p1;
  input >> p1;
  if (!input)
  {
    std::cerr << "invalid OFF file" << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  // read in target model, apply transformation and compare
  input.open(fname2);
  if (!input)
  {
    std::cerr << "unable to open file" << std::endl;
    return EXIT_FAILURE;
  }
  Polyhedron p2;
  input >> p2;
  CGAL::Polygon_mesh_processing::transform(res2, p2);
  Mesh m;
  CGAL::copy_face_graph(p2, m);

  std::cout << "Construct AABB tree...";
  Tree tree(faces(p1).first, faces(p1).second, p1);
  std::cout << "done." << std::endl;

  Mesh::Property_map<vertex_descriptor, FT> distance;
  bool created;
  boost::tie(distance, created) = m.add_property_map<vertex_descriptor, FT>("v:quality", 0.0);
  assert(created);
  std::cout << "computing distances for the tooth..." << std::endl;
  for (vertex_descriptor vi : m.vertices())
  {
    if (vertexSet2.find(vi) != vertexSet2.end())
    {
      // if on the tooth, compute the distance
      Point_3 point = m.point(vi);
      // perform distance query
      FT dist = CGAL::sqrt(tree.squared_distance(point));
      if (dist < 0.2)
      { // tolerance
        dist = 0;
      }
      // get sign through ray casting (random vector)
      typedef Tree::size_type size_type;
      Vector_3 random_vec = random_vector_3();
      Ray_3 ray(point, random_vec);
      size_type nbi = tree.number_of_intersected_primitives(ray);
      FT sign((nbi & 1) ? 1 : -1);
      dist = sign * dist;

      distance[vi] = dist;
      // std::cout << dist << std::endl;
    }
    else
    {
      // if not on the tooth, set value to 0
      distance[vi] = 0;
    }
  }
  std::ofstream out("diff.ply");
  bool wrote = CGAL::write_ply(out, m);
  assert(wrote);

  std::cout << "done." << std::endl;
}
