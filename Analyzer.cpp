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
#include <numeric>
#include <stdlib.h>
#include <time.h>

#include <CGAL/intersections.h>

#include "Analyzer.h"
#include "utility.h"
#include "objects.h"
#include "result.h"

// additional headers for analyze
#include <CGAL/squared_distance_3.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/property_map.h>
#include <CGAL/enum.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/read_off_points.h>

#include <CGAL/Polygon_mesh_processing/measure.h> // face_area()

#include <iterator>
#include <string>

#include "Roughness.h"

#include <QFuture> // for async computation
#include <QtConcurrent/QtConcurrent> // for async computation
#include <QFile>

using namespace mycode;

int Analyzer::analyze()
{
  emit updateProgressBar(10);
  init();
  emit updateProgressBar(20);
  compute_shoulder_width();
  emit updateProgressBar(33);
  compute_axial_wall_height();
  emit updateProgressBar(46);
  compute_taper();
  emit updateProgressBar(59);
  compute_margin_depth();
  emit updateProgressBar(72);
  compute_occlusal_reduction();
  emit updateProgressBar(85);
  compute_gingival_extension();
  emit updateProgressBar(90);
  // compute_roughness();
  emit updateProgressBar(100);
  return 0;
}

void Analyzer::init()
{
  std::ifstream input(this->param.studentModel);
  CGAL::read_off(input, this->student_model);
  input.close();

  input.open(this->param.studentModel);
  input >> this->student_model_poly;
  input.close();

  input.open(this->param.originalModel);
  CGAL::read_off(input, this->original_model);
  input.close();

  input.open(this->param.originalModel);
  input >> this->original_model_poly;
  input.close();

  /* compute alignment matrix
   * and align student model with original model */
  K::Aff_transformation_3 transformation = this->compute_alignment_matrix();
  if (debug) {
    std::cout << "transformation" << transformation << std::endl;
  }
  CGAL::Polygon_mesh_processing::transform(transformation, this->student_model);
  CGAL::Polygon_mesh_processing::transform(transformation, this->student_model_poly);

  if (debug) {
    // write transformed model to output
    std::string output_name = "aligned1.ply";
    std::ofstream out(output_name);
    CGAL::write_ply(out, student_model);
  }

  student_tree.rebuild(faces(this->student_model_poly).first, faces(this->student_model_poly).second, this->student_model_poly);
  original_tree.rebuild(faces(this->original_model_poly).first, faces(this->original_model_poly).second, this->original_model_poly);

  if (param.divisionEnabled) {
    std::vector<mycode::Point_3> temp;
    readpp(temp, param.studentCenterPoint);
    if (temp.size() != 1) {
      std::cerr << "student model center point is not a single point, size: " << temp.size() << std::endl;
      return;
    }
    student_model_center = temp[0].transform(transformation);
    temp.clear();
    readpp(temp, param.studentMidpoint);
    if (temp.size() != 1) {
      std::cerr << "student model mid point is not a single point, size: " << temp.size() << std::endl;
      return;
    }
    student_model_midpoint = temp[0].transform(transformation);
 }

  /* read sets of points on student model from file
   * and transform these points according to the alignment matrix */
  readpp(margin_points, param.studentMarginPoints);
  for (auto& point : margin_points) {
    point = point.transform(transformation);
  }
  readpp(axial_points, param.studentAxialPoints);
  for (auto& point : axial_points) {
    point = point.transform(transformation);
  }
  readpp(occlusal_points, param.studentOcclusalPoints);
  for (auto& point : occlusal_points) {
    point = point.transform(transformation);
  }
  readpp(gingiva_points, param.studentGingivaPoints);
  for (auto& point : gingiva_points) {
    point = point.transform(transformation);
  }
}

K::Aff_transformation_3 Analyzer::compute_alignment_matrix()
{
  if (param.useManualTransform) {
    return K::Aff_transformation_3(param.transformMatrix[0][0],param.transformMatrix[0][1],param.transformMatrix[0][2],param.transformMatrix[0][3],
        param.transformMatrix[1][0],param.transformMatrix[1][1],param.transformMatrix[1][2],param.transformMatrix[1][3],
        param.transformMatrix[2][0],param.transformMatrix[2][1],param.transformMatrix[2][2],param.transformMatrix[2][3]);
  } else {
    std::vector<Pwn> pwns1, pwns2;
    std::ifstream input(param.studentModel);
    CGAL::read_off_points(input, std::back_inserter(pwns1),
                CGAL::parameters::point_map (Point_map()).
                normal_map (Normal_map()));
    input.close();

    input.open(param.originalModel);
    CGAL::read_off_points(input, std::back_inserter(pwns2),
                CGAL::parameters::point_map (Point_map()).
                normal_map (Normal_map()));
    input.close();
    K::Aff_transformation_3 res1 =
        std::get<0>(
            CGAL::OpenGR::compute_registration_transformation(pwns1, pwns2,
                                                              params::point_map(Point_map()).normal_map(Normal_map()),
                                                              params::point_map(Point_map()).normal_map(Normal_map())));

    K::Aff_transformation_3 res2 =
        std::get<0>(
            CGAL::pointmatcher::compute_registration_transformation(pwns1, pwns2,
                                                                    params::point_map(Point_map()).normal_map(Normal_map()),
                                                                    params::point_map(Point_map()).normal_map(Normal_map()).transformation(res1)));
    return res2;
  }
}

// constructing lines from points
void Analyzer::construct_lines(std::vector<mycode::Point_3> &points, std::vector<mycode::Segment_2> &lines)
{
  auto start = points.begin();
  mycode::Point_2 start_2(start->x(), start->z());
  for (auto p = points.begin(); p != points.end(); p++)
  {
    mycode::Point_2 p_2(p->x(), p->z());
    auto q = std::next(p, 1);
    if (q == points.end()) {
      lines.push_back(mycode::Segment_2(p_2, start_2));
    } else {
      mycode::Point_2 q_2(q->x(), q->z());
      lines.push_back(mycode::Segment_2(p_2, q_2));
    }
  }
}

// return true if a point is within a closed 2D segment set
bool Analyzer::within_lines(mycode::Point_2 &point, std::vector<mycode::Segment_2> &lines)
{
  mycode::Vector_2 random_vec = random_vector_2();
  mycode::Ray_2 ray(point, random_vec);
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
  // used when division enabled
  std::vector<mycode::FT> tapers_lingual;
  std::vector<mycode::FT> tapers_buccal;
  std::vector<mycode::FT> tapers_mesial;
  std::vector<mycode::FT> tapers_distal;

  // used when division disabled
  std::vector<mycode::FT> tapers;

  /* compute toc
   * first we compute the "average" point of occlusal points */
  int count = 0;
  mycode::FT x_avg = 0;
  mycode::FT z_avg = 0;
  for (auto op : occlusal_points)
  {
    x_avg += op.x();
    z_avg += op.z();
    count++;
  }
  std::cerr << "count = " << count << std::endl;
  x_avg /= count;
  z_avg /= count;

  /* for every occlusal point, we calculate its corresponding toc */
  for (auto target_point : occlusal_points)
  {
    /* compute far point */
    mycode::FT max_angle = -1;
    mycode::Point_3 far_point;
    for (auto op : occlusal_points)
    {
      mycode::FT angle = CGAL::approximate_angle(mycode::Point_3(target_point.x(), 0, target_point.z()),
                                                 mycode::Point_3(x_avg, 0, z_avg),
                                                 mycode::Point_3(op.x(), 0, op.z()));
      if (angle > max_angle)
      {
        far_point = op;
        max_angle = angle;
      }
    }

    /* compute axial wall approx of target point */
    mycode::Point_3 vertical_point;
    mycode::FT min_angle = 180;
    for (auto ap : axial_points)
    {
      mycode::FT angle = CGAL::approximate_angle(mycode::Vector_3(ap, target_point), mycode::Vector_3(0, 1, 0));
      if (angle < min_angle)
      {
        vertical_point = ap;
        min_angle = angle;
      }
    }
    /* use mid point for approximating axial wall slope */
    mycode::Point_3 mid_point = student_tree.closest_point(CGAL::midpoint(target_point, vertical_point));
    mycode::Vector_3 v1 = mycode::Vector_3(vertical_point, mid_point);

    /* compute axial wall approx of far point */
    min_angle = 180;
    for (auto ap : axial_points)
    {
      mycode::FT angle = CGAL::approximate_angle(mycode::Vector_3(ap, far_point), mycode::Vector_3(0, 1, 0));
      if (angle < min_angle)
      {
        vertical_point = ap;
        min_angle = angle;
      }
    }

    /* use midpoint for approximating axial wall slope */
    mid_point = student_tree.closest_point(CGAL::midpoint(far_point, vertical_point));
    mycode::Vector_3 v2 = mycode::Vector_3(vertical_point, mid_point);

    /* toc is the angle formed by the axial wall approx of target point and far point */
    mycode::FT toc = CGAL::approximate_angle(v1, v2);
    mycode::FT taper = toc/2;

    if (param.divisionEnabled) {
      std::cerr << "checking region" << std::endl;
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
  }

  /* report stats to the console */
  emit msgToConsole("=============");
  emit msgToConsole("=== TAPER ===");
  emit msgToConsole("=============");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    report_stats(&student_result.taper[0], tapers_lingual);
    emit msgToConsole("=== Buccal ===");
    report_stats(&student_result.taper[1], tapers_buccal);
    emit msgToConsole("=== Mesial ===");
    report_stats(&student_result.taper[2], tapers_mesial);
    emit msgToConsole("=== Distal ===");
    report_stats(&student_result.taper[3], tapers_distal);
  } else {
    report_stats(&student_result.taper[0], tapers);
  }
}

void Analyzer::compute_occlusal_reduction()
{
  std::vector<mycode::FT> occlusal_reductions_lingual;
  std::vector<mycode::FT> occlusal_reductions_buccal;
  std::vector<mycode::FT> occlusal_reductions_mesial;
  std::vector<mycode::FT> occlusal_reductions_distal;

  std::vector<mycode::FT> occlusal_reductions;

  std::unordered_set<mycode::vertex_descriptor> points_on_occlusal;
  select_occlusal_points(points_on_occlusal);
  for (auto& vi : points_on_occlusal) {
    mycode::Point_3 p = student_model.point(vi);
    mycode::FT dist = CGAL::sqrt(original_tree.squared_distance(p));
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
  }

  /* report stats to the console */
  emit msgToConsole("==========================");
  emit msgToConsole("=== OCCLUSAL REDUCTION ===");
  emit msgToConsole("==========================");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    report_stats(&student_result.occlusal_reduction[0], occlusal_reductions_lingual);
    emit msgToConsole("=== Buccal ===");
    report_stats(&student_result.occlusal_reduction[1], occlusal_reductions_buccal);
    emit msgToConsole("=== Mesial ===");
    report_stats(&student_result.occlusal_reduction[2], occlusal_reductions_mesial);
    emit msgToConsole("=== Distal ===");
    report_stats(&student_result.occlusal_reduction[3], occlusal_reductions_distal);
  } else {
    report_stats(&student_result.occlusal_reduction[0], occlusal_reductions);
  }
}

void Analyzer::compute_margin_depth()
{
  std::vector<mycode::FT> margin_depths_lingual;
  std::vector<mycode::FT> margin_depths_buccal;
  std::vector<mycode::FT> margin_depths_mesial;
  std::vector<mycode::FT> margin_depths_distal;

  std::vector<mycode::FT> margin_depths;

  for (auto& p : margin_points) {
    mycode::Point_3 point_half_mm_above(p.x(), p.y() + 0.5, p.z());
    mycode::Point_3 student_point = student_tree.closest_point(point_half_mm_above);
    mycode::Point_3 original_point = original_tree.closest_point(point_half_mm_above);
    mycode::FT dist = CGAL::sqrt(CGAL::squared_distance(student_point, original_point));
    if (param.divisionEnabled) {
      switch (region_of(p)) {
        case Lingual:
          margin_depths_lingual.push_back(dist);
          break;
        case Buccal:
          margin_depths_buccal.push_back(dist);
          break;
        case Mesial:
          margin_depths_mesial.push_back(dist);
          break;
        case Distal:
          margin_depths_distal.push_back(dist);
          break;
      }
    } else {
      margin_depths.push_back(dist);
    }
  }

  /* report stats to the console */
  emit msgToConsole("====================");
  emit msgToConsole("=== MARGIN DEPTH ===");
  emit msgToConsole("====================");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    report_stats(&student_result.margin_depth[0], margin_depths_lingual);
    emit msgToConsole("=== Buccal ===");
    report_stats(&student_result.margin_depth[1], margin_depths_buccal);
    emit msgToConsole("=== Mesial ===");
    report_stats(&student_result.margin_depth[2], margin_depths_mesial);
    emit msgToConsole("=== Distal ===");
    report_stats(&student_result.margin_depth[3], margin_depths_distal);
  } else {
    report_stats(&student_result.margin_depth[0], margin_depths);
  }
}

void Analyzer::compute_gingival_extension()
{
  std::vector<mycode::FT> y_mp_lingual;
  std::vector<mycode::FT> y_mp_buccal;
  std::vector<mycode::FT> y_mp_mesial;
  std::vector<mycode::FT> y_mp_distal;

  std::vector<mycode::FT> y_gp_lingual;
  std::vector<mycode::FT> y_gp_buccal;
  std::vector<mycode::FT> y_gp_mesial;
  std::vector<mycode::FT> y_gp_distal;

  std::vector<mycode::FT> y_mp;
  std::vector<mycode::FT> y_gp;

  for (auto& mp : margin_points) {
    if (param.divisionEnabled) {
      switch (region_of(mp)) {
        case Lingual:
          y_mp_lingual.push_back(mp.y());
          break;
        case Buccal:
          y_mp_buccal.push_back(mp.y());
          break;
        case Mesial:
          y_mp_mesial.push_back(mp.y());
          break;
        case Distal:
          y_mp_distal.push_back(mp.y());
          break;
      }
    } else {
      y_mp.push_back(mp.y());
    }
  }

  for (auto& gp : gingiva_points) {
    if (param.divisionEnabled) {
      switch (region_of(gp)) {
        case Lingual:
          y_gp_lingual.push_back(gp.y());
          break;
        case Buccal:
          y_gp_buccal.push_back(gp.y());
          break;
        case Mesial:
          y_gp_mesial.push_back(gp.y());
          break;
        case Distal:
          y_gp_distal.push_back(gp.y());
          break;
      }
    } else {
      y_gp.push_back(gp.y());
    }
  }

  /* report stats to the console */

  mycode::FT average_y_mp, average_y_gp, average_gingival_extension, max_gingival_extension, min_gingival_extension;

  emit msgToConsole("==========================");
  emit msgToConsole("=== GINGIVAL EXTENSION ===");
  emit msgToConsole("==========================");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    average_y_mp = accumulate(y_mp_lingual.begin(), y_mp_lingual.end(), 0.0)/y_mp_lingual.size();
    average_y_gp = accumulate(y_gp_lingual.begin(), y_gp_lingual.end(), 0.0)/y_gp_lingual.size();
    average_gingival_extension = average_y_mp - average_y_gp;
//    max_gingival_extension = *max_element(y_mp_lingual.begin(), y_mp_lingual.end()) - *min_element(y_gp_lingual.begin(), y_gp_lingual.end());
//    min_gingival_extension = *min_element(y_mp_lingual.begin(), y_mp_lingual.end()) - *max_element(y_gp_lingual.begin(), y_gp_lingual.end());
    emit msgToConsole(QString("average = %1").arg(average_gingival_extension));
//    emit msgToConsole(QString("max = %1").arg(max_gingival_extension));
//    emit msgToConsole(QString("min = %1").arg(min_gingival_extension));

    emit msgToConsole("=== Buccal ===");
    average_y_mp = accumulate(y_mp_buccal.begin(), y_mp_buccal.end(), 0.0)/y_mp_buccal.size();
    average_y_gp = accumulate(y_gp_buccal.begin(), y_gp_buccal.end(), 0.0)/y_gp_buccal.size();
    average_gingival_extension = average_y_mp - average_y_gp;
//    max_gingival_extension = *max_element(y_mp_buccal.begin(), y_mp_buccal.end()) - *min_element(y_gp_buccal.begin(), y_gp_buccal.end());
//    min_gingival_extension = *min_element(y_mp_buccal.begin(), y_mp_buccal.end()) - *max_element(y_gp_buccal.begin(), y_gp_buccal.end());
    emit msgToConsole(QString("average = %1").arg(average_gingival_extension));
//    emit msgToConsole(QString("max = %1").arg(max_gingival_extension));
//    emit msgToConsole(QString("min = %1").arg(min_gingival_extension));

    emit msgToConsole("=== Mesial ===");
    average_y_mp = accumulate(y_mp_mesial.begin(), y_mp_mesial.end(), 0.0)/y_mp_mesial.size();
    average_y_gp = accumulate(y_gp_mesial.begin(), y_gp_mesial.end(), 0.0)/y_gp_mesial.size();
    average_gingival_extension = average_y_mp - average_y_gp;
//    max_gingival_extension = *max_element(y_mp_mesial.begin(), y_mp_mesial.end()) - *min_element(y_gp_mesial.begin(), y_gp_mesial.end());
//    min_gingival_extension = *min_element(y_mp_mesial.begin(), y_mp_mesial.end()) - *max_element(y_gp_mesial.begin(), y_gp_mesial.end());
    emit msgToConsole(QString("average = %1").arg(average_gingival_extension));
//    emit msgToConsole(QString("max = %1").arg(max_gingival_extension));
//    emit msgToConsole(QString("min = %1").arg(min_gingival_extension));

    emit msgToConsole("=== Distal ===");
    average_y_mp = accumulate(y_mp_distal.begin(), y_mp_distal.end(), 0.0)/y_mp_distal.size();
    average_y_gp = accumulate(y_gp_distal.begin(), y_gp_distal.end(), 0.0)/y_gp_distal.size();
    average_gingival_extension = average_y_mp - average_y_gp;
//    max_gingival_extension = *max_element(y_mp_distal.begin(), y_mp_distal.end()) - *min_element(y_gp_distal.begin(), y_gp_distal.end());
//    min_gingival_extension = *min_element(y_mp_distal.begin(), y_mp_distal.end()) - *max_element(y_gp_distal.begin(), y_gp_distal.end());
    emit msgToConsole(QString("average = %1").arg(average_gingival_extension));
//    emit msgToConsole(QString("max = %1").arg(max_gingival_extension));
//    emit msgToConsole(QString("min = %1").arg(min_gingival_extension));
  } else {
    average_y_mp = accumulate(y_mp.begin(), y_mp.end(), 0.0)/y_mp.size();
    average_y_gp = accumulate(y_gp.begin(), y_gp.end(), 0.0)/y_gp.size();
    average_gingival_extension = average_y_mp - average_y_gp;
//    max_gingival_extension = *max_element(y_mp.begin(), y_mp.end()) - *min_element(y_gp.begin(), y_gp.end());
//    min_gingival_extension = *min_element(y_mp.begin(), y_mp.end()) - *max_element(y_gp.begin(), y_gp.end());
    emit msgToConsole(QString("average = %1").arg(average_gingival_extension));
//    emit msgToConsole(QString("max = %1").arg(max_gingival_extension));
//    emit msgToConsole(QString("min = %1").arg(min_gingival_extension));
  }
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
  }

  /* report stats to the console */
  emit msgToConsole("======================");
  emit msgToConsole("=== SHOULDER WIDTH ===");
  emit msgToConsole("======================");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    report_stats(&student_result.shoulder_width[0], shoulder_widths_lingual);
    emit msgToConsole("=== Buccal ===");
    report_stats(&student_result.shoulder_width[1], shoulder_widths_buccal);
    emit msgToConsole("=== Mesial ===");
    report_stats(&student_result.shoulder_width[2], shoulder_widths_mesial);
    emit msgToConsole("=== Distal ===");
    report_stats(&student_result.shoulder_width[3], shoulder_widths_distal);
  } else {
    report_stats(&student_result.shoulder_width[0], shoulder_widths);
  }
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
  }

  /* report stats to the console */
  emit msgToConsole("=========================");
  emit msgToConsole("=== AXIAL WALL HEIGHT ===");
  emit msgToConsole("=========================");
  if (param.divisionEnabled) {
    emit msgToConsole("=== Lingual ===");
    report_stats(&student_result.axial_wall_height[0], axial_wall_height_lingual);
    emit msgToConsole("=== Buccal ===");
    report_stats(&student_result.axial_wall_height[1], axial_wall_height_buccal);
    emit msgToConsole("=== Mesial ===");
    report_stats(&student_result.axial_wall_height[2], axial_wall_height_mesial);
    emit msgToConsole("=== Distal ===");
    report_stats(&student_result.axial_wall_height[3], axial_wall_height_distal);
  } else {
    report_stats(&student_result.axial_wall_height[0], axial_wall_height);
  }
}

void Analyzer::compute_roughness()
{
  emit msgToConsole("computing roughness for the shoulder, this might take about 1 minute...");

  std::unordered_set<mycode::vertex_descriptor> pointsOnShoulder;
  std::unordered_set<mycode::vertex_descriptor> pointsOnAxialWall;
  select_shoulder_points(pointsOnShoulder);
  select_axial_wall_points(pointsOnAxialWall);

  PolyhedronPtr poly = PolyhedronPtr(new Polyhedron_enriched()); // it was called Polyhedron, Polyhedron_enriched is just an alias to resolve typedef issues
  poly->load_mesh_off(param.studentModel);
  mycode::Mesh m1; // m1 will store roughness on shoulder
  mycode::Mesh m2; // m2 will store roughness on axial wall

  std::ifstream input(this->param.studentModel);
  CGAL::read_off(input, m1);
  input.close();

  input.open(this->param.studentModel);
  CGAL::read_off(input, m2);
  input.close();

  poly->Normalise();
  poly->compute_bounding_box();
  poly->compute_normals();
  poly->compute_type();
  poly->calc_nb_components();
  poly->calc_nb_boundaries();

  double epsilon = 0.01;

  CRoughness<Polyhedron_enriched> roughness(poly.get());
  roughness.compute_Roughness(2 * epsilon, epsilon);

  mycode::Mesh::Property_map<mycode::vertex_descriptor, double> roughness_on_shoulder;
  bool created;
  boost::tie(roughness_on_shoulder, created) = m1.add_property_map<mycode::vertex_descriptor, double>("v:quality", 0.0);
  assert(created);
  mycode::Mesh::Property_map<mycode::vertex_descriptor, double> roughness_on_axial_wall;
  boost::tie(roughness_on_axial_wall, created) = m2.add_property_map<mycode::vertex_descriptor, double>("v:quality", 0.0);
  assert(created);

  // roughness of points in four regions
  std::vector<mycode::FT> shoulder_roughness_lingual;
  std::vector<mycode::FT> shoulder_roughness_buccal;
  std::vector<mycode::FT> shoulder_roughness_mesial;
  std::vector<mycode::FT> shoulder_roughness_distal;

  std::vector<mycode::FT> axial_wall_roughness_lingual;
  std::vector<mycode::FT> axial_wall_roughness_buccal;
  std::vector<mycode::FT> axial_wall_roughness_mesial;
  std::vector<mycode::FT> axial_wall_roughness_distal;

  // average roughness in four sectors
  double avg_roughness_shoulder[4] = {};
  double avg_roughness_axial_wall[4] = {};

  mycode::vertex_descriptor vi = *(m1.vertices().begin());
  for (Vertex_iterator pVertex = poly->vertices_begin(); pVertex != poly->vertices_end(); pVertex++)
  {
    if (pointsOnShoulder.find(vi) != pointsOnShoulder.end()) {
      // add roughness quality on the shoulder
      mycode::Point_3 p = this->student_model.point(vi);
      roughness_on_shoulder[vi] = pVertex->Roughness();
      switch (region_of(p)) {
        case Lingual:
          shoulder_roughness_lingual.push_back(pVertex->Roughness());
          break;
        case Buccal:
          shoulder_roughness_buccal.push_back(pVertex->Roughness());
          break;
        case Mesial:
          shoulder_roughness_mesial.push_back(pVertex->Roughness());
          break;
        case Distal:
          shoulder_roughness_distal.push_back(pVertex->Roughness());
          break;
      }
    } else if (pointsOnAxialWall.find(vi) != pointsOnAxialWall.end()) {
      // add roughness quality on the axial wall
      mycode::Point_3 p = this->student_model.point(vi);
      roughness_on_axial_wall[vi] = pVertex->Roughness();
      switch (region_of(p)) {
        case Lingual:
          axial_wall_roughness_lingual.push_back(pVertex->Roughness());
          break;
        case Buccal:
          axial_wall_roughness_buccal.push_back(pVertex->Roughness());
          break;
        case Mesial:
          axial_wall_roughness_mesial.push_back(pVertex->Roughness());
          break;
        case Distal:
          axial_wall_roughness_distal.push_back(pVertex->Roughness());
          break;
      }
    }
    vi++;
  }
  avg_roughness_shoulder[0] = std::accumulate(shoulder_roughness_lingual.begin(), shoulder_roughness_lingual.end(), 0.0) / shoulder_roughness_lingual.size();
  avg_roughness_shoulder[1] = std::accumulate(shoulder_roughness_buccal.begin(), shoulder_roughness_buccal.end(), 0.0) / shoulder_roughness_buccal.size();
  avg_roughness_shoulder[2] = std::accumulate(shoulder_roughness_mesial.begin(), shoulder_roughness_mesial.end(), 0.0) / shoulder_roughness_mesial.size();
  avg_roughness_shoulder[3] = std::accumulate(shoulder_roughness_distal.begin(), shoulder_roughness_distal.end(), 0.0) / shoulder_roughness_distal.size();

  avg_roughness_axial_wall[0] = std::accumulate(axial_wall_roughness_lingual.begin(), axial_wall_roughness_lingual.end(), 0.0) / axial_wall_roughness_lingual.size();
  avg_roughness_axial_wall[1] = std::accumulate(axial_wall_roughness_buccal.begin(), axial_wall_roughness_buccal.end(), 0.0) / axial_wall_roughness_buccal.size();
  avg_roughness_axial_wall[2] = std::accumulate(axial_wall_roughness_mesial.begin(), axial_wall_roughness_mesial.end(), 0.0) / axial_wall_roughness_mesial.size();
  avg_roughness_axial_wall[3] = std::accumulate(axial_wall_roughness_distal.begin(), axial_wall_roughness_distal.end(), 0.0) / axial_wall_roughness_distal.size();

  std::ofstream out("roughness_shoulder.ply");
  bool wrote = CGAL::write_ply(out, m1);
  if (!wrote)
  {
    std::cerr << "write failed" << std::endl;
  }

  std::ofstream out2("roughness_axial_wall.ply");
  wrote = CGAL::write_ply(out2, m2);
  if (!wrote)
  {
    std::cerr << "write failed" << std::endl;
  }

  /* report stats to the console */
  emit msgToConsole("=============================");
  emit msgToConsole("=== ROUGHNESS ON SHOULDER ===");
  emit msgToConsole("=============================");
  emit msgToConsole("=== Lingual ===");
  report_stats(&student_result.roughness_shoulder[0], shoulder_roughness_lingual);
  emit msgToConsole("=== Buccal ===");
  report_stats(&student_result.roughness_shoulder[1], shoulder_roughness_buccal);
  emit msgToConsole("=== Mesial ===");
  report_stats(&student_result.roughness_shoulder[2], shoulder_roughness_mesial);
  emit msgToConsole("=== Distal ===");
  report_stats(&student_result.roughness_shoulder[3], shoulder_roughness_distal);

  emit msgToConsole("===============================");
  emit msgToConsole("=== ROUGHNESS ON AXIAL WALL ===");
  emit msgToConsole("===============================");
  emit msgToConsole("=== Lingual ===");
  report_stats(&student_result.roughness_axial_wall[0], axial_wall_roughness_lingual);
  emit msgToConsole("=== Buccal ===");
  report_stats(&student_result.roughness_axial_wall[1], axial_wall_roughness_buccal);
  emit msgToConsole("=== Mesial ===");
  report_stats(&student_result.roughness_axial_wall[2], axial_wall_roughness_mesial);
  emit msgToConsole("=== Distal ===");
  report_stats(&student_result.roughness_axial_wall[3], axial_wall_roughness_distal);
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

void Analyzer::report_stats(double* avg, const std::vector<mycode::FT> &values)
{
  mycode::FT average = accumulate(values.begin(), values.end(), 0.0)/values.size();
  *avg = average; /* save the average value into student_result */
  emit msgToConsole(QString("number of values in region = %1").arg(values.size()));
  emit msgToConsole(QString("average = %1").arg(average));
  emit msgToConsole(QString("max = %1").arg(*max_element(values.begin(), values.end())));
  emit msgToConsole(QString("min = %1").arg(*min_element(values.begin(), values.end())));
}

void Analyzer::feedback()
{
  return;
}
