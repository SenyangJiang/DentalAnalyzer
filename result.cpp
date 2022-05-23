#include "result.h"

QString Stats::to_csv() {
  return QString("%1, %2, %3,\n").arg(this->max).arg(this->min).arg(this->avg);
}

QString points_to_csv(std::vector<std::pair<mycode::Point_3, double>>& data) {
  QString res = "";
  for (auto &pair : data) {
    res.append(QString("\"(%1,%2,%3)\",").arg(pair.first.x()).arg(pair.first.y()).arg(pair.first.z()));
  }
  res.append("\n");
  return res;
}

QString values_to_csv(std::vector<std::pair<mycode::Point_3, double>>& data) {
  QString res = "";
  for (auto &pair : data) {
    res.append(QString("%1,").arg(pair.second));
  }
  res.append("\n");
  return res;
}
