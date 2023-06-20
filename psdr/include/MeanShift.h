#pragma once 
//#include "defs.h"
#include "defs_cgal_ui.h"
class MeanShift {
public:


	MeanShift() { }

    std::pair<double, double> meanshift(const std::pair<double, double> &point_original, const std::vector<std::pair<double, double>> &all_point, const std::vector<double> &points_number, const std::vector<int> &if_consider,
		double kernel_bandwidth,
		double EPSILON);

	double meanshift_oned(const double &distance_original, const std::vector<double> &all_distances, const std::vector<double> &points_number, const std::vector<int> &if_consider,
		double kernel_bandwidth,
		double EPSILON);

private:
	double euclidean_distance(const std::pair<double, double> &point_a, const std::pair<double, double> &point_b);
	double euclidean_distance_sqr(const std::pair<double, double> &point_a, const std::pair<double, double> &point_b);

	std::pair<double, double> shift_point(const std::pair<double, double> &point,
		const std::vector<std::pair<double, double>> &points, const std::vector<double> &points_number,
		double kernel_bandwidth, const std::vector<int> &if_consider);

	double euclidean_distance_oned(const double &point_a, const double &point_b);
	double euclidean_distance_sqr_oned(const double &point_a, const double &point_b);
    double shift_distance(const double &distance_input,
		const std::vector<double> &distance_inputs, const std::vector<double> &points_number,
		double kernel_bandwidth, const std::vector<int> &if_consider);

};
