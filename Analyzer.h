#ifndef ANALYZER_H_INCLUDED
#define ANALYZER_H_INCLUDED

#include <QObject>

#include <vector>

#include "objects.h"
#include "result.h"
#include "parameter.h"

class Analyzer : public QObject
{
    Q_OBJECT

public:
    // constructor
    Analyzer(Parameter param);

    // compute function called by mainwindow
    int analyze();

    Parameter param;

    bool status_done;

    Result student_result;

signals:
    void msgToConsole(QString msg);

    void updateProgressBar(int val);

    void alertToWindow(QString msg);

private:
    // initialize private variables
    bool init();

    // parse .off file to mesh and polyhedral surface
    bool read_off(string file, mycode::Mesh &mesh, mycode::Polyhedron &poly);

    // read in points in .pp file
    bool readpp(std::vector<mycode::Point_3> &points, std::string filename);

    // constructing 2d lines from points
    void construct_lines(std::vector<mycode::Point_3> &points, std::vector<mycode::Segment_2> &lines);

    // return true if a 2D point is within a 2D segment set
    bool within_lines(mycode::Point_2 &p, std::vector<mycode::Segment_2> &lines);

    // compute taper (half of toc)
    void compute_taper();

    // compute occlusal reduction
    void compute_occlusal_reduction();

    // compute margin depth
    void compute_margin_depth();

    // compute gingival extension
    void compute_gingival_extension();

    // computes shoulder width
    void compute_shoulder_width();

    // compute axial wall height
    void compute_axial_wall_height();

    // determine what region the point belongs to
    Region region_of(mycode::Point_3 point);

    // select vertices on the occlusal
    void select_occlusal_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet);

    // select vertices on the shoulder
    void select_shoulder_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet);

    // select vertices on the axial wall
    void select_axial_wall_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet);

    // select vertices on tooth
    void select_tooth_points(std::unordered_set<mycode::vertex_descriptor> &vertexSet);

    // write stats to console
    void report_stats(QString metric, Stats* stats_arr,
                      const std::vector<mycode::FT> &values,
                      const std::vector<mycode::FT> &values_lingual,
                      const std::vector<mycode::FT> &values_buccal,
                      const std::vector<mycode::FT> &values_mesial,
                      const std::vector<mycode::FT> &values_distal);

    // report avg, min, max of an array of numbers, and record them in stats struct
    void compute_stats(Stats &stats, const std::vector<mycode::FT> &values);

    bool debug = false; /* debug flag for testing */

    // private variables to read in
    mycode::Mesh student_model;
    mycode::Mesh original_model;
    mycode::Polyhedron student_model_poly;
    mycode::Polyhedron original_model_poly;
    mycode::Tree student_tree;
    mycode::Tree original_tree;
    mycode::Point_3 student_model_center;
    mycode::Point_3 student_model_midpoint;
    std::vector<mycode::Point_3> occlusal_points;
    std::vector<mycode::Point_3> axial_points;
    std::vector<mycode::Point_3> margin_points;
    std::vector<mycode::Point_3> gingiva_points;
};
#endif
