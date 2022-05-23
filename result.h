#ifndef RESULT_H
#define RESULT_H

#include <vector>
#include <utility>
#include <objects.h>
#include <QString>

class Stats {
public:
  double max;
  double min;
  double avg;

  QString to_csv();
};

struct Result {
  Stats shoulder_width_stats[4]; /* Lingual, Buccal, Mesial, Distal */
  Stats taper_stats[4];
  Stats axial_wall_height_stats[4];
  Stats margin_depth_stats[4];
  Stats occlusal_reduction_stats[4];
  Stats gingival_extension_stats[4];

  std::vector<std::pair<mycode::Point_3, double>> occlusal_reduction_data;
  std::vector<std::pair<mycode::Point_3, double>> taper_data;
  std::vector<std::pair<mycode::Point_3, double>> margin_depth_data;
  std::vector<std::pair<mycode::Point_3, double>> gingival_extension_data;
  std::vector<std::pair<mycode::Point_3, double>> shoulder_width_data;
  std::vector<std::pair<mycode::Point_3, double>> axial_wall_height_data;
};

QString points_to_csv(std::vector<std::pair<mycode::Point_3, double>>& data);

QString values_to_csv(std::vector<std::pair<mycode::Point_3, double>>& data);

#endif // RESULT_H
