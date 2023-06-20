#include <math.h>
#include "MeanShift.h"

double MeanShift::euclidean_distance(const std::pair<double, double> &point_a, const std::pair<double, double> &point_b) {
	double total = (point_a.first - point_b.first)*(point_a.first - point_b.first) + (point_a.second - point_b.second)*(point_a.second - point_b.second);
	return std::sqrt(total);
}

double MeanShift::euclidean_distance_sqr(const std::pair<double, double> &point_a, const std::pair<double, double> &point_b) {
	double total = (point_a.first - point_b.first)*(point_a.first - point_b.first) + (point_a.second - point_b.second)*(point_a.second - point_b.second);
	return total;
}

double gaussian_kernel(double distance, double kernel_bandwidth) {
	double temp = exp(-1.0 / 2.0 * (distance *distance) / (kernel_bandwidth*kernel_bandwidth));
	return temp;
}



std::pair<double, double> MeanShift::shift_point(const std::pair<double, double> &point,
	const std::vector<std::pair<double, double>> &points, const std::vector<double> &points_number,
	double kernel_bandwidth, const std::vector<int> &if_consider) {

	std::pair<double, double> shift_vector = std::make_pair(0, 0);
	double total_weight = 0;
	for (int i = 0; i < points.size(); i++) {
		if (if_consider[i] != -1) continue;
		std::pair<double, double> temp_point = points[i];
		double distance = euclidean_distance(point, temp_point);

		double gauss_weight = gaussian_kernel(distance, kernel_bandwidth);
		shift_vector.first += temp_point.first*gauss_weight*points_number[i];
		shift_vector.second += temp_point.second*gauss_weight*points_number[i];



		total_weight += gauss_weight * points_number[i];
	}
	shift_vector.first = shift_vector.first / total_weight;
	shift_vector.second = shift_vector.second / total_weight;

	return shift_vector;

}

std::pair<double, double> MeanShift::meanshift(const std::pair<double, double> &point_original, const std::vector<std::pair<double, double>> &all_point, const std::vector<double> &points_number, const std::vector<int> &if_consider,
	double kernel_bandwidth,
	double EPSILON) {
	const double EPSILON_SQR = EPSILON * EPSILON;

	std::pair<double, double> shifted_points = point_original;

	bool keep_moving = true;

	do {

		keep_moving = false;
		std::pair<double, double> new_point = shift_point(shifted_points, all_point, points_number, kernel_bandwidth, if_consider);
		double shift_distance_sqr = euclidean_distance_sqr(new_point, shifted_points);

		if (shift_distance_sqr >= EPSILON_SQR) {
			shifted_points = new_point;
			keep_moving = true;
		}

	} while (keep_moving);
	return shifted_points;
}


double MeanShift::meanshift_oned(const double &distance_original, const std::vector<double> &all_distances, const std::vector<double> &points_number, const std::vector<int> &if_consider,
	double kernel_bandwidth,
	double EPSILON) {
	const double EPSILON_SQR = EPSILON * EPSILON;

	double shifted_distnaces = distance_original;

	bool keep_moving = true;

	do {

		keep_moving = false;
		double new_distance = shift_distance(shifted_distnaces, all_distances, points_number, kernel_bandwidth, if_consider);
		double shift_distance_sqr = euclidean_distance_sqr_oned(new_distance, shifted_distnaces);

		if (shift_distance_sqr >= EPSILON_SQR) {
			shifted_distnaces = new_distance;
			keep_moving = true;
		}

	} while (keep_moving);
	return shifted_distnaces;
}


double MeanShift::shift_distance(const double &distance_input,
	const std::vector<double> &distance_inputs, const std::vector<double> &points_number,
	double kernel_bandwidth, const std::vector<int> &if_consider) {

	double shift_vector = 0;
	double total_weight = 0;
	for (int i = 0; i < distance_inputs.size(); i++) {
		if (if_consider[i] != -1) continue;
		double temp_point = distance_inputs[i];
		double distance = euclidean_distance_oned(distance_input, temp_point);

		double gauss_weight = gaussian_kernel(distance, kernel_bandwidth);
		shift_vector += temp_point * gauss_weight*points_number[i];




		total_weight += gauss_weight * points_number[i];
	}
	shift_vector = shift_vector / total_weight;


	return shift_vector;

}

double MeanShift::euclidean_distance_oned(const double &point_a, const double &point_b) {
	double total = (point_a - point_b)*(point_a - point_b);
	return std::sqrt(total);
}


double MeanShift::euclidean_distance_sqr_oned(const double &point_a, const double &point_b) {
	double total = (point_a - point_b)*(point_a - point_b);
	return total;
}
