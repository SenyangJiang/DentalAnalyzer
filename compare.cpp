#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/property_map.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/pointmatcher/register_point_sets.h>
#include <CGAL/pointmatcher/compute_registration_transformation.h>
#include <CGAL/OpenGR/compute_registration_transformation.h>

#include <CGAL/Polygon_mesh_processing/transform.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <stdlib.h>
#include <time.h>

// for registration
typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::FT FT;
typedef K::Ray_3 Ray_3;
typedef std::pair<Point_3, Vector_3> Pwn;
typedef CGAL::First_of_pair_property_map<Pwn> Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn> Normal_map;

// for constructing AABB tree and distance computation
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace params = CGAL::parameters;

// utilities for generating random vector
FT random_in(const float a, const float b)
{
    float r = rand() / (float)RAND_MAX;
    return (FT)(a + (b - a) * r);
}


Vector_3 random_vector()
{
    FT x = random_in(0.0,1.0);
    FT y = random_in(0.0,1.0);
    FT z = random_in(0.0,1.0);
    return Vector_3(x,y,z);
}

// main function
int main(int argc, const char** argv) {
  const char* fname1 = "data/reference.off";
  const char* fname2 = "data/1.off";
  std::vector<Pwn> pwns1, pwns2;
  std::ifstream input(fname1);
  if (!input ||
      !CGAL::read_off_points(input, std::back_inserter(pwns1),
			     CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
			     normal_map(Normal_map())))
    {
      std::cerr << "Error: cannot read file " << fname1 << std::endl;
      return EXIT_FAILURE;
    }
  input.close();
  input.open(fname2);
  if (!input ||
      !CGAL::read_off_points(input, std::back_inserter(pwns2),
			     CGAL::parameters::point_map(Point_map()).
			     normal_map(Normal_map())))
    {
      std::cerr << "Error: cannot read file " << fname2 << std::endl;
      return EXIT_FAILURE;
    }
  input.close();
    
  std::cerr << "Computing registration transformation using OpenGR Super4PCS.." << std::endl;
  // First, compute registration transformation using OpenGR Super4PCS
  K::Aff_transformation_3 res1 =
    std::get<0>( // get first of pair, which is the transformation
		CGAL::OpenGR::compute_registration_transformation
		(pwns1, pwns2,
		 params::point_map(Point_map()).normal_map(Normal_map()),
		 params::point_map(Point_map()).normal_map(Normal_map()))
		 );
  std::cerr << "Computing registration transformation using PointMatcher ICP, "
          << "taking transformation computed by OpenGR Super4PCS as initial transformation.." << std::endl;
  K::Aff_transformation_3 res2 =
    std::get<0>(
		CGAL::pointmatcher::compute_registration_transformation
		(pwns1, pwns2,
		 params::point_map(Point_map()).normal_map(Normal_map()),
		 params::point_map(Point_map()).normal_map(Normal_map()).
		 transformation(res1) /* initial transform for pwns2.
				       * a transform returned from a coarse registration algorithm.
				       * */
		 ));
  
  input.open(fname1);
  if (!input) {
    std::cerr << "unable to open file" << std::endl;
    return EXIT_FAILURE;
  }
  Polyhedron p1;
  input >> p1;
  if (!input) {
    std::cerr << "invalid OFF file" << std::endl;
    return EXIT_FAILURE;
  }
  input.close();
  input.open(fname2);
  if (!input) {
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

  Mesh::Property_map<vertex_descriptor,FT> distance;
  bool created;
  boost::tie(distance, created) = m.add_property_map<vertex_descriptor,FT>("v:quality",0.0);
  assert(created);
  std::cout << "computing distances..." << std::endl;
  for (vertex_descriptor vi : m.vertices()) {
    Point_3 point = m.point(vi);
    // perform distance query
    FT dist = CGAL::sqrt(tree.squared_distance(point));
    // get sign through ray casting (random vector)
    srand (time(NULL));
    typedef Tree::size_type size_type;
    Vector_3 random_vec = random_vector();
    Ray_3 ray(point, random_vec);
    size_type nbi = tree.number_of_intersected_primitives(ray);
    FT sign ( (nbi&1) == 0 ? 1 : -1);
    dist = sign * dist;

    distance[vi] = dist;
    // std::cout << dist << std::endl;
  }
  std::ofstream out("diff.ply");
  bool wrote = CGAL::write_ply(out, m);
  assert(wrote);
  
  std::cout << "done." << std::endl;
}
