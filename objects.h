#ifndef OBJECTS_H_INCLUDED
#define OBJECTS_H_INCLUDED

/*
    This file includes all typedefs used in our project
*/

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/property_map.h>

enum Region {
  Lingual,
  Buccal,
  Mesial,
  Distal
};

namespace mycode
{
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point_3;
    typedef K::Segment_3 Segment_3;
    typedef K::Line_3 Line_3;
    typedef K::Vector_3 Vector_3;
    typedef K::FT FT;
    typedef K::Ray_3 Ray_3;

    typedef K::Ray_2 Ray_2;
    typedef K::Vector_2 Vector_2;
    typedef K::Segment_2 Segment_2;
    typedef K::Point_2 Point_2;

    typedef CGAL::Surface_mesh<Point_3> Mesh;
    typedef Mesh::Vertex_index vertex_descriptor;
    typedef Mesh::Face_index face_descriptor;
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    typedef std::pair<Point_3, Vector_3> Pwn;
    typedef CGAL::First_of_pair_property_map<Pwn> Point_map;
    typedef CGAL::Second_of_pair_property_map<Pwn> Normal_map;

    namespace params = CGAL::parameters;
} // namespace mycode

#endif
