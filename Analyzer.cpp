#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/property_map.h>
#include <CGAL/enum.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <numeric>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iterator>
#include <string>

#include "Analyzer.h"
#include "utility.h"
#include "objects.h"
#include "result.h"

#include <QFuture> // for async computation
#include <QtConcurrent/QtConcurrent> // for async computation
#include <QFile>
using namespace mycode;

// constructor
Analyzer::Analyzer(Parameter param)
{
  status_done = false;
  this->param = param;
}

// compute function called by mainwindow
int Analyzer::analyze()
{
  emit updateProgressBar(0);
  if (!init()) {
    status_done = true;
    return 1;
  }
  emit updateProgressBar(20);
  compute_shoulder_width();
  emit updateProgressBar(33);
  compute_axial_wall_height();
  emit updateProgressBar(46);
  compute_taper();
  //emit updateProgressBar(59);
  //compute_margin_depth();
  emit updateProgressBar(72);
  compute_occlusal_reduction();
  emit updateProgressBar(85);
  compute_gingival_extension();
  emit updateProgressBar(100);

  /* mark status as done */
  status_done = true;
  return 0;
}

// initialize private variables
bool Analyzer::init()
{
  /* read models */
  if (!read_off(this->param.studentModel, this->student_model, this->student_model_poly)) {
    return false;
  }

  if (!read_off(this->param.originalModel, this->original_model, this->original_model_poly)) {
    return false;
  }

  /* build AABB tree for distance query */
  student_tree.rebuild(faces(this->student_model_poly).first, faces(this->student_model_poly).second, this->student_model_poly);
  original_tree.rebuild(faces(this->original_model_poly).first, faces(this->original_model_poly).second, this->original_model_poly);

  /* read optional center & mid point */
  if (param.divisionEnabled) {
    std::vector<Point_3> temp;
    if (!readpp(temp, param.studentCenterPoint)) {
      return false;
    }
    if (temp.size() != 1) {
      emit alertToWindow(QString("student model center point is not a single point, size: %1").arg(temp.size()));
      return false;
    }
    temp.clear();
    if (!readpp(temp, param.studentMidpoint)) {
      return false;
    }
    if (temp.size() != 1) {
      emit alertToWindow(QString("student model mid point is not a single point, size: %1").arg(temp.size()));
      return false;
    }
 }

  /* read sets of points on student model from file */
  readpp(margin_points, param.studentMarginPoints);
  if (margin_points.size() < 20) {
    emit alertToWindow(QString("please specify more then 20 points for margin line, current size: %1").arg(margin_points.size()));
    return false;
  }
  readpp(axial_points, param.studentAxialPoints);
  if (axial_points.size() < 20) {
    emit alertToWindow(QString("please specify more then 20 points for axial line, current size: %1").arg(axial_points.size()));
    return false;
  }
  readpp(occlusal_points, param.studentOcclusalPoints);
  if (occlusal_points.size() < 20) {
    emit alertToWindow(QString("please specify more then 20 points for occlusal line, current size: %1").arg(occlusal_points.size()));
    return false;
  }
  readpp(gingiva_points, param.studentGingivaPoints);
  if (gingiva_points.size() < 20) {
    emit alertToWindow(QString("please specify more then 20 points for gingiva line, current size: %1").arg(gingiva_points.size()));
    return false;
  }

  return true;
}

// parse .off file to mesh and polyhedral surface
bool Analyzer::read_off(string file, Mesh &mesh, Polyhedron &poly) {
  std::ifstream input;
  input.open(file);
  if (!input.good()) {
    emit alertToWindow(QString("model file '%1' cannot be read").arg(QString::fromStdString(file)));
    return false;
  }
  if (!CGAL::read_off(input, mesh)) {
    emit alertToWindow(QString("cannot parse model file '%1' to mesh").arg(QString::fromStdString(file)));
    return false;
  }
  input.close();

  input.open(file);
  if (!input.good()) {
    emit alertToWindow(QString("model file '%1' cannot be read").arg(QString::fromStdString(file)));
    return false;
  }
  if (!CGAL::read_off(input, poly)) {
    emit alertToWindow(QString("cannot parse model file '%1' to polyhedral surface").arg(QString::fromStdString(file)));
    return false;
  }
  input.close();
  return true;
}

// function for reading in points
bool Analyzer::readpp(std::vector<Point_3> &points, std::string filename)
{
  std::ifstream input(filename);
  if (!input.good())
  {
    emit alertToWindow(QString("point file '%1' cannot be read").arg(QString::fromStdString(filename)));
    return false;
  }
  double x;
  double y;
  double z;
  std::string line;
  std::size_t start;
  std::size_t end;
  while (std::getline(input, line))
  {
    end = line.find("<point");
    if (end != std::string::npos)
    {
      while (1)
      {
        start = line.find("\"", end + 1);
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
      points.push_back(Point_3(x, y, z));
    }
  }
  input.close();
  return true;
}

// constructing 2d lines from points
void Analyzer::construct_lines(std::vector<Point_3> &points, std::vector<Segment_2> &lines)
{
  auto start = points.begin();
  Point_2 start_2(start->x(), start->z());
  for (auto p = points.begin(); p != points.end(); p++)
  {
    Point_2 p_2(p->x(), p->z());
    auto q = std::next(p, 1);
    if (q == points.end()) {
      lines.push_back(Segment_2(p_2, start_2));
    } else {
      Point_2 q_2(q->x(), q->z());
      lines.push_back(Segment_2(p_2, q_2));
    }
  }
}

// return true if a point is within a closed 2D segment set
bool Analyzer::within_lines(Point_2 &point, std::vector<Segment_2> &lines)
{
  Vector_2 random_vec = random_vector_2();
  Ray_2 ray(point, random_vec);
  int count = 0;
  for (auto l : lines)
  {
    if (intersection(ray, l))
      count++;
  }
  return count % 2 == 1;
}

// compute taper (half of toc)
void Analyzer::compute_taper()
{
  /* used when division enabled */
  std::vector<double> tapers_lingual;
  std::vector<double> tapers_buccal;
  std::vector<double> tapers_mesial;
  std::vector<double> tapers_distal;

  /* used when division disabled */
  std::vector<double> tapers;

  /* first we compute the "average" point of occlusal points
   * which will be the center of the circle */
  int count = 0;
  double x_avg = 0;
  double z_avg = 0;
  for (auto op : occlusal_points)
  {
    x_avg += op.x();
    z_avg += op.z();
    count++;
  }
  x_avg /= count;
  z_avg /= count;

  /* for every occlusal point, we calculate its corresponding toc */
  for (auto target_point : occlusal_points)
  {
    /* compute far point */
    double max_angle = -1;
    Point_3 far_point;
    for (auto op : occlusal_points)
    {
      double angle = CGAL::approximate_angle(Point_3(target_point.x(), 0, target_point.z()),
                                                 Point_3(x_avg, 0, z_avg),
                                                 Point_3(op.x(), 0, op.z()));
      if (angle > max_angle)
      {
        far_point = op;
        max_angle = angle;
      }
    }

    /* compute axial wall approx of target point */
    Point_3 vertical_point;
    double min_angle = 180;
    for (auto ap : axial_points)
    {
      double angle = CGAL::approximate_angle(Vector_3(ap, target_point), Vector_3(0, 1, 0));
      if (angle < min_angle)
      {
        vertical_point = ap;
        min_angle = angle;
      }
    }
    /* use mid point for approximating axial wall slope */
    Point_3 mid_point = student_tree.closest_point(CGAL::midpoint(target_point, vertical_point));
    Vector_3 v1(vertical_point, mid_point);

    /* similarly, compute axial wall approx of far point */
    min_angle = 180;
    for (auto ap : axial_points)
    {
      double angle = CGAL::approximate_angle(Vector_3(ap, far_point), Vector_3(0, 1, 0));
      if (angle < min_angle)
      {
        vertical_point = ap;
        min_angle = angle;
      }
    }

    /* similarly, use midpoint for approximating axial wall slope */
    mid_point = student_tree.closest_point(CGAL::midpoint(far_point, vertical_point));
    mycode::Vector_3 v2 = mycode::Vector_3(vertical_point, mid_point);

    /* toc is the angle formed by the axial wall approx of target point and far point */
    mycode::FT toc = CGAL::approximate_angle(v1, v2);
    mycode::FT taper = toc/2;

    if (param.divisionEnabled) {
      switch (region_of(target_point)) {
        case Lingual:
          tapers_lingual.push_back(taper);
          break;
        case Buccal:
          tapers_buccal.push_back(taper);
          break;
        case Mesial:
          tapers_mesial.push_back(taper);
          break;
        case Distal:
          tapers_distal.push_back(taper);
          break;
      }
    } else {
      tapers.push_back(taper);
    }

    student_result.taper_data.push_back(std::make_pair(target_point, taper));
  }

  /* report stats to the console */
  report_stats("TAPER", student_result.taper_stats, tapers, tapers_lingual, tapers_buccal, tapers_mesial, tapers_distal);
}

void Analyzer::compute_occlusal_reduction()
{
  std::vector<mycode::FT> occlusal_reductions_lingual;
  std::vector<mycode::FT> occlusal_reductions_buccal;
  std::vector<mycode::FT> occlusal_reductions_mesial;
  std::vector<mycode::FT> occlusal_reductions_distal;

  std::vector<mycode::FT> occlusal_reductions;

  std::unordered_set<mycode::vertex_descriptor> points_on_occlusal;
  std::unordered_set<mycode::vertex_descriptor> points_on_original;
   /* new update:   */
  occlusal_projected_points(points_on_original);
  for (auto& vi : points_on_original) {
    mycode::Point_3 p = original_model.point(vi);
    mycode::FT dist = CGAL::sqrt(student_tree.squared_distance(p));

//  select_occlusal_points(points_on_occlusal);
//  for (auto& vi : points_on_occlusal) {
//    mycode::Point_3 p = student_model.point(vi);
//    mycode::FT dist = CGAL::sqrt(original_tree.squared_distance(p));
    if (param.divisionEnabled) {
      switch (region_of(p)) {
        case Lingual:
          occlusal_reductions_lingual.push_back(dist);
          break;
        case Buccal:
          occlusal_reductions_buccal.push_back(dist);
          break;
        case Mesial:
          occlusal_reductions_mesial.push_back(dist);
          break;
        case Distal:
          occlusal_reductions_distal.push_back(dist);
          break;
      }
    } else {
      occlusal_reductions.push_back(dist);
    }

    student_result.occlusal_reduction_data.push_back(std::make_pair(p, dist));
  }

  /* report stats to the console */
   report_stats("OCCLUSAL REDUCTION", student_result.occlusal_reduction_stats, occlusal_reductions, occlusal_reductions_lingual, occlusal_reductions_buccal, occlusal_reductions_mesial, occlusal_reductions_distal);
}

//void Analyzer::compute_margin_depth()
//{
//  std::vector<mycode::FT> margin_depths_lingual;
//  std::vector<mycode::FT> margin_depths_buccal;
//  std::vector<mycode::FT> margin_depths_mesial;
//  std::vector<mycode::FT> margin_depths_distal;

//  std::vector<mycode::FT> margin_depths;

//  for (auto p : margin_points) {
//    mycode::Point_3 point_half_mm_above(p.x(), p.y() + 0.5, p.z());
//    mycode::Point_3 student_point = student_tree.closest_point(point_half_mm_above);
//    mycode::Point_3 original_point = original_tree.closest_point(point_half_mm_above);
//    mycode::FT dist = CGAL::sqrt(CGAL::squared_distance(student_point, original_point));
//    if (param.divisionEnabled) {
//      switch (region_of(p)) {
//        case Lingual:
//          margin_depths_lingual.push_back(dist);
//          break;
//        case Buccal:
//          margin_depths_buccal.push_back(dist);
//          break;
//        case Mesial:
//          margin_depths_mesial.push_back(dist);
//          break;
//        case Distal:
//          margin_depths_distal.push_back(dist);
//          break;
//      }
//    } else {
//      margin_depths.push_back(dist);
//    }

//    student_result.margin_depth_data.push_back(std::make_pair(p, dist));
//  }

//  report_stats("MARGIN DEPTH", student_result.margin_depth_stats, margin_depths, margin_depths_lingual, margin_depths_buccal, margin_depths_mesial,margin_depths_distal);
//}

void Analyzer::compute_gingival_extension()
{
  std::vector<mycode::FT> gingival_extensions_lingual;
  std::vector<mycode::FT> gingival_extensions_buccal;
  std::vector<mycode::FT> gingival_extensions_mesial;
  std::vector<mycode::FT> gingival_extensions_distal;

  std::vector<mycode::FT> gingival_extensions;

  for (auto mp : margin_points)
  {
    mycode::FT min_angle = 180;
    mycode::Point_3 vertical_point;
    for (auto gp : gingiva_points)
    {
      mycode::Vector_3 v(mp, gp);
      mycode::FT angle = std::min(CGAL::approximate_angle(v, mycode::Vector_3(0, 1, 0)), CGAL::approximate_angle(v, mycode::Vector_3(0, -1, 0)));
      if (angle < min_angle)
      {
        vertical_point = gp;
        min_angle = angle;
      }
    }
    mycode::FT height_diff = mp.y() - vertical_point.y();
    if (param.divisionEnabled) {
      switch (region_of(mp)) {
        case Lingual:
          gingival_extensions_lingual.push_back(height_diff);
          break;
        case Buccal:
          gingival_extensions_buccal.push_back(height_diff);
          break;
        case Mesial:
          gingival_extensions_mesial.push_back(height_diff);
          break;
        case Distal:
          gingival_extensions_distal.push_back(height_diff);
          break;
      }
    } else {
      gingival_extensions.push_back(height_diff);
    }

    student_result.gingival_extension_data.push_back(std::make_pair(mp, height_diff));
  }

  /* report stats to the console */
  report_stats("GINGIVAL EXTENSION", student_result.gingival_extension_stats, gingival_extensions, gingival_extensions_lingual, gingival_extensions_buccal, gingival_extensions_mesial, gingival_extensions_distal);
}

void Analyzer::compute_shoulder_width()
{
  std::vector<mycode::FT> shoulder_widths_lingual;
  std::vector<mycode::FT> shoulder_widths_buccal;
  std::vector<mycode::FT> shoulder_widths_mesial;
  std::vector<mycode::FT> shoulder_widths_distal;

  std::vector<mycode::FT> shoulder_widths;

  for (auto ap : axial_points)
  {
    mycode::FT min_dist = -1;
    for (auto mp : margin_points)
    {
      mycode::FT dist = CGAL::sqrt(CGAL::squared_distance(ap, mp));
      if (min_dist < 0 || dist < min_dist)
      {
        min_dist = dist;
      }
    }
    if (param.divisionEnabled) {
      switch (region_of(ap)) {
        case Lingual:
          shoulder_widths_lingual.push_back(min_dist);
          break;
        case Buccal:
          shoulder_widths_buccal.push_back(min_dist);
          break;
        case Mesial:
          shoulder_widths_mesial.push_back(min_dist);
          break;
        case Distal:
          shoulder_widths_distal.push_back(min_dist);
          break;
      }
    } else {
      shoulder_widths.push_back(min_dist);
    }

    student_result.shoulder_width_data.push_back(std::make_pair(ap, min_dist));
  }

  /* report stats to the console */
  report_stats("SHOULDER WIDTH", student_result.shoulder_width_stats, shoulder_widths, shoulder_widths_lingual, shoulder_widths_buccal, shoulder_widths_mesial, shoulder_widths_distal);
}

void Analyzer::compute_axial_wall_height()
{
  std::vector<mycode::FT> axial_wall_height_lingual;
  std::vector<mycode::FT> axial_wall_height_buccal;
  std::vector<mycode::FT> axial_wall_height_mesial;
  std::vector<mycode::FT> axial_wall_height_distal;

  std::vector<mycode::FT> axial_wall_height;

  for (auto op : occlusal_points)
  {
    mycode::FT min_angle = 180;
    mycode::Point_3 vertical_point;
    for (auto ap : axial_points)
    {
      mycode::FT angle = CGAL::approximate_angle(mycode::Vector_3(ap, op), mycode::Vector_3(0, 1, 0));
      if (angle < min_angle)
      {
        vertical_point = ap;
        min_angle = angle;
      }
    }
    mycode::FT height = CGAL::sqrt(CGAL::squared_distance(op, vertical_point));
    if (param.divisionEnabled) {
      switch (region_of(op)) {
        case Lingual:
          axial_wall_height_lingual.push_back(height);
          break;
        case Buccal:
          axial_wall_height_buccal.push_back(height);
          break;
        case Mesial:
          axial_wall_height_mesial.push_back(height);
          break;
        case Distal:
          axial_wall_height_distal.push_back(height);
          break;
      }
    } else {
      axial_wall_height.push_back(height);
    }

    student_result.axial_wall_height_data.push_back(std::make_pair(op, height));
  }

  /* report stats to the console */
  report_stats("AXIAL WALL HEIGHT", student_result.axial_wall_height_stats, axial_wall_height, axial_wall_height_lingual, axial_wall_height_buccal, axial_wall_height_mesial, axial_wall_height_distal);
}

Region Analyzer::region_of(mycode::Point_3 point)
{
  /* find the center of tooth */
  double min_x = margin_points.begin()->x();
  double max_x = margin_points.begin()->x();
  double min_z = margin_points.begin()->z();
  double max_z = margin_points.begin()->z();
  for (auto& p : margin_points)
  {
    if (p.x() < min_x)
    {
      min_x = p.x();
    }
    else if (p.x() > max_x)
    {
      max_x = p.x();
    }

    if (p.z() < min_z)
    {
      min_z = p.z();
    }
    else if (p.z() > max_z)
    {
      max_z = p.z();
    }
  }

  double center_x = (min_x + max_x) / 2;
  double center_z = (min_z + max_z) / 2;
  mycode::Point_2 tooth_center(center_x, center_z);

  /* find the closest point to model_center as the center of Lingual Surface */
  mycode::Point_2 model_center(student_model_center.x(), student_model_center.z());
  mycode::Point_2 model_midpoint(student_model_midpoint.x(), student_model_midpoint.z());
  mycode::Point_2 lingual_center;
  double min_dist = 100;
  for (auto& point : margin_points) {
    mycode::Point_2 p(point.x(), point.z());
    double dist = CGAL::sqrt(CGAL::squared_distance(model_center, p));
    if (dist < min_dist) {
      lingual_center = p;
      min_dist = dist;
    }
  }

  double reference_angle = CGAL::approximate_angle(mycode::Point_3(model_midpoint.x(), 0, model_midpoint.y()),
                                         mycode::Point_3(model_center.x(), 0, model_center.y()),
                                         mycode::Point_3(lingual_center.x(), 0, lingual_center.y()));

  double compare_angle = CGAL::approximate_angle(mycode::Point_3(model_midpoint.x(), 0, model_midpoint.y()),
                                         mycode::Point_3(model_center.x(), 0, model_center.y()),
                                         mycode::Point_3(point.x(), 0, point.z()));
  double angle_with_lingual_center = CGAL::approximate_angle(mycode::Point_3(lingual_center.x(), 0, lingual_center.y()),
                                         mycode::Point_3(tooth_center.x(), 0, tooth_center.y()),
                                         mycode::Point_3(point.x(), 0, point.z()));
  if (angle_with_lingual_center < 45) {
    return Lingual;
  } else if (angle_with_lingual_center < 135 && compare_angle < reference_angle) {
    return Mesial;
  } else if (angle_with_lingual_center < 135 && compare_angle > reference_angle) {
    return Distal;
  } else {
    return Buccal;
  }
}


void Analyzer::occlusal_projected_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet){  // upward projection from occlusal points to obtain points on original

    std::unordered_set<mycode::vertex_descriptor> points_on_occlusal;
    select_occlusal_points(points_on_occlusal);

     for (auto& p : points_on_occlusal){
        mycode::Point_3 occ_p = student_model.point(p);

       for(mycode::vertex_descriptor ori : original_model.vertices())
       {
          mycode::Point_3 ori_p = original_model.point(ori);

          if( std::abs(ori_p.x() - occ_p.x() ) < 0.09 && std::abs(ori_p.z() - occ_p.z()) < 0.09 && (ori_p.y() > occ_p.y()))        // *new update* select vertex points from the orginal model if x and z are similar with the occlusal points
          {
              vertexSet.insert(ori);
          }
       }
    }

}
void Analyzer::select_occlusal_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet)
{
  /* calculate y average of occlusal points */
  mycode::FT occlusal_avg_y = 0;
  int count = 0;


  for (auto p = occlusal_points.begin(); p != occlusal_points.end(); p++)
  {
    occlusal_avg_y += p->y();
    count++;
  }
  occlusal_avg_y /= count;

  std::vector<mycode::Segment_2> occlusal_lines;
  construct_lines(occlusal_points, occlusal_lines);

  for (mycode::vertex_descriptor vi : student_model.vertices())
  {
    mycode::Point_3 p = student_model.point(vi);
    mycode::Point_2 point(p.x(), p.z());

    if (within_lines(point, occlusal_lines) && (std::abs(p.y() - occlusal_avg_y) < 2)) /* This '2' is pretty random, not sure how to select points with similar height as the occlusal line */
    {
      vertexSet.insert(vi);
    }
  }

}

void Analyzer::select_shoulder_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet)
{
  /* calculate y average of margin points */
  mycode::FT margin_avg_y = 0;
  int count = 0;
  for (auto p = margin_points.begin(); p != margin_points.end(); p++)
  {
    margin_avg_y += p->y();
    count++;
  }
  margin_avg_y /= count;

  std::vector<mycode::Segment_2> margin_lines;
  std::vector<mycode::Segment_2> axial_lines;
  construct_lines(margin_points, margin_lines);
  construct_lines(axial_points, axial_lines);

  for (mycode::vertex_descriptor vi : student_model.vertices())
  {
    mycode::Point_3 p = student_model.point(vi);
    mycode::Point_2 point(p.x(), p.z());

    /* the point is on the shoulder if
     * 1. it is inside the margin lines
     * 2. it is outside the axial lines */
    if (within_lines(point, margin_lines) && !within_lines(point, axial_lines) && (std::abs(p.y() - margin_avg_y) < 1))
    {
      vertexSet.insert(vi);
    }
  }
}

void Analyzer::select_axial_wall_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet)
{
    std::unordered_set<mycode::vertex_descriptor> points_on_tooth;
    std::unordered_set<mycode::vertex_descriptor> points_on_occlusal;
    std::unordered_set<mycode::vertex_descriptor> points_on_shoulder;
    select_tooth_points(points_on_tooth);
    select_occlusal_points(points_on_occlusal);
    select_shoulder_points(points_on_shoulder);

    for (vertex_descriptor vi : points_on_tooth) {
      if ((points_on_occlusal.find(vi) == points_on_occlusal.end()) && (points_on_shoulder.find(vi) == points_on_shoulder.end())) {
        vertexSet.insert(vi);
      }
    }
}

void Analyzer::select_tooth_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet)
{
  /* calculate the mean value of y */
  mycode::FT avg_y = 0;
  int count = 0;
  for (auto p = margin_points.begin(); p != margin_points.end(); p++)
  {
    avg_y += p->y();
    count++;
  }
  avg_y /= count;
  std::vector<mycode::Segment_2> margin_lines;
  construct_lines(margin_points, margin_lines);

  for (vertex_descriptor vi : student_model.vertices())
  {
    mycode::Point_3 p = student_model.point(vi);
    mycode::Point_2 point(p.x(), p.z());

    if (within_lines(point, margin_lines) && (std::abs(p.y() - avg_y) < 10))
    {
      vertexSet.insert(vi);
    }
  }
}

void Analyzer::report_stats(QString metric, Stats* stats_arr,
                            const std::vector<mycode::FT> &values,
                            const std::vector<mycode::FT> &values_lingual,
                            const std::vector<mycode::FT> &values_buccal,
                            const std::vector<mycode::FT> &values_mesial,
                            const std::vector<mycode::FT> &values_distal)
{
  emit msgToConsole(QString("<h2 style=\"color:#FF0c32;\">%1</h2>").arg(metric));
  if (param.divisionEnabled) {
    emit msgToConsole("<h3>Lingual</h3>");
    compute_stats(stats_arr[0], values_lingual);
    emit msgToConsole("<h3>Buccal</h3>");
    compute_stats(stats_arr[1], values_buccal);
    emit msgToConsole("<h3>Mesial</h3>");
    compute_stats(stats_arr[2], values_mesial);
    emit msgToConsole("<h3>Distal</h3>");
    compute_stats(stats_arr[3], values_distal);
  } else {
    compute_stats(stats_arr[0], values);
  }
}

void Analyzer::compute_stats(Stats &stats, const std::vector<mycode::FT> &values) {
  double average = accumulate(values.begin(), values.end(), 0.0)/values.size();
  double max = *max_element(values.begin(), values.end());
  double min = *min_element(values.begin(), values.end());
  stats.avg = average;
  stats.max = max;
  stats.min = min;
  emit msgToConsole(QString("number of values in region = <b>%1</b>").arg(values.size()));
  emit msgToConsole(QString("average = <b>%1</b>").arg(average));
  emit msgToConsole(QString("max = <b>%1</b>").arg(max));
  emit msgToConsole(QString("min = <b>%1</b>").arg(min));
}
