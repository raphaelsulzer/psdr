#pragma once

#include <fstream>

#include "defs.h"
#include "defs_cgal_ui.h"
#include "MeanShift.h"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Shape_detection.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

using namespace std;


typedef pair<Inexact_Point_3, Inexact_Vector_3> Point_with_normal;
typedef vector<Point_with_normal> Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef CGAL::Shape_detection::Shape_detection_traits
<CGAL::Exact_predicates_inexact_constructions_kernel, Pwn_vector, Point_map, Normal_map> Traits;
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


class Shape_Detector
{

	struct Add_Comparator {
        bool operator()(pair<int, double> p1, pair<int, double> p2) {

			return p1.second > p2.second;
		}
	};
	struct Add_Comparator_normal {
        bool operator()(pair<int, double> p1, pair<int, double> p2) {

			return p1.second < p2.second;
		}
	};

	struct Remove_Comparator {
        bool operator()(pair<int, double> p1, pair<int, double> p2) {

			return p1.second < p2.second;
		}
	};
	struct Remove_Comparator_normal {
        bool operator()(pair<int, double> p1, pair<int, double> p2) {

			return p1.second > p2.second;
		}
	};

	struct Weight_Comparator_with_energy {
        bool operator()(pair<pair<int, vector<int>>, vector<double>> p1, pair<pair<int, vector<int>>, vector<double>> p2) {

			return p1.second[0] > p2.second[0];
		}
	};

public:
    Shape_Detector();

    std::shared_ptr<spdlog::logger> _logger;

    int load_points(const string _filename);

protected:
	bool load_ply();

	bool load_vg();
    bool load_npz();


protected:
	void set_extrema();

public:

    int set_detection_parameters(
            int _min_points,
            double _epsilon,
            int _knn,
            double _normal_threshold);

    void set_regularization_parameters(
		double _lambda,
		double _tolerance_angle,
		double _tolerance_coplanarity);

    void set_knn(int _knn);

    void set_discretization_parameters(double _discretization_angle,
                                       double _discretization_distance_ratio);

	void detect_shapes();

	void load_shapes();

	void regularize_shapes();

    void refine_shapes(int mi);

	void planar_shape_detection_L1();
	void get_last_information();
	void test_connected_primitives();
	void get_distance_diviation();
	void get_good_points_normal();

protected:
	void compute_average_spacing_and_k_nearest_neighbors();

	void reshuffle_points_based_on_planarity();

	void do_region_growing();

public:
    void export_region_growing_results(const string & filename);

    void set_path_clusters(const string & filename);

	void load_region_growing_results();

protected:
	void detect_planes(bool read_clusters);

	void regularize_planes();

	void discretize_planes();

    bool inliers_arent_aligned(const vector<int> & indices);

	void get_coverage_and_mean_error();
    double add_changed_error_normal(int id, vector<int> point_ids);
	double energy_changed_second_normal(double dis, double numb, int id_plane);
	void get_bad_points_normal();
	void local_operators_normal();
    bool convert_form_merge_normal(int i, int j, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_split_normal(int i, vector<vector<int>>& max_list_vector, vector<vector<int>>& min_list_vector, pair<pair<int, vector<int>>, vector<double>> & one_element);
    bool convert_form_exclude_normal(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_insert_normal(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool separate_two_out_normal(int id, vector<int> & max_list, vector<int> & min_list, double & dif);
	double energy_changed_normal(double dis, double numb, int id_plane);
    bool merge_distance_changed_normal_with_epsilon(int i, int j, double & dif, vector<int> & move_ids);
	double energy_changed_normal_merge(double dis, double numb, int nmoves, int id_plane_1, int id_plane_2);
    double remove_changed_error_normal(int id, vector<int> point_ids);
    bool update_bad_points_normal(int id_shape, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool update_good_points_normal(int id_shape, pair<pair<int, vector<int>>, vector<double>>& one_element);
	void get_distance_diviation_show_merge_info(double t);
	void transfer_operator_normal();
	bool test_if_connected(int i, int j);
	void get_distance_diviation_show_normal_info(double t);
	void get_coverage_and_mean_error_pure();
	void planar_shape_detection_l2();
	void get_good_points();
    double add_changed_error(int id, vector<int> point_ids);
	double energy_changed_second(double dis, double numb, int id_plane);
	int change_of_degrees_of_freedom_after_add_remove(int id_plane);
	int change_of_degrees_of_freedom_after_mergy(int id_1, int id_2);
	int change_of_degrees_of_freedom_after_split(int id_plane);
	void get_bad_points();
    double remove_changed_error(int id, vector<int> point_ids);
	void local_operators();
    void update_regular_relations_after_merge(int id1, int id2, vector<int> respective_planes);
	void update_regular_done_group_by_planes();
    bool update_good_points(int id_shape, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool update_bad_points(int id_shape, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_insert_l2(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_regularization_l2(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_exclude_l2(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
    bool convert_form_split_l2(int i, vector<vector<int>>& max_list_vector, vector<vector<int>>& min_list_vector, pair<pair<int, vector<int>>, vector<double>> & one_element);
    bool merge_distance_changed_with_epsilon(int i, int j, double & dif, vector<int> & move_ids);
    bool convert_form_merge_l2(int i, int j, pair<pair<int, vector<int>>, vector<double>>& one_element);
	double energy_changed_merge(double dis, double numb, int nmoves, int id_plane_1, int id_plane_2);
    bool separate_two_out(int id, vector<int> & max_list, vector<int> & min_list, double & dif);
	double energy_changed(double dis, double numb, int id_plane);
	void transfer_operator_l2();
	int get_the_recent_degrees_of_freedom();
	void planar_shape_detection_hybrid();
	void transfer_operator_normal_for_hybrid();
    void calculate_center();
	void show_result(double t);
    pair<double, double> get_gauss_sphere_value(const Inexact_Vector_3 & v);
public:

//    void set_default_parameters();

	void set_primitives();

	void test_regularity_planes();

	void clear_primitives();

	void set_stop_iteration(bool cc);

	void set_constraint(bool cc);

	void set_lambda_r(double db);

	void set_lambda_c(double db);

	void set_lambda_regularity(double db);

	void set_lambda_fidelity(double db);

	void set_weight_m(int wm);

    void set_metric_type(int m);

    void set_max_iter(int mi);

	int get_number_of_non_coplanar_planes();

	std::tuple<Inexact_Point_3, Inexact_Vector_3, CGAL::Color, CGAL::Color> get_vertex(int point_index);

    void get_inliers(size_t primitive_index, vector<Point_with_normal> & inliers);
    void get_inliers(size_t primitive_index, vector<size_t> & inliers);

	size_t get_total_number_of_inliers();

	int get_cloud_size() const;

	void discretize_vector(const double ds, const Inexact_Vector_3 & v, int & k_longitude, int & k_latitude);

    void make_3d_histogram(std::map<pair<int, int>, list<int> > & atlas, const double ds, vector<pair<Inexact_Vector_3, double> > & H);
		
    void print_3d_histogram(const string & filename, vector<pair<Inexact_Vector_3, double> > & H);

    void to_vg(const string& directory, string & basename, int X, int Y, int Z);

    void to_vg(const string& filename);

    void to_npz(const string& filename);



	Pwn_vector::const_iterator points_begin() const;

	Pwn_vector::const_iterator points_end() const;
	void plane_parellel_group_mean_shift();


	void initialize_parallel_clusters();
	void initialize_orthogonal_clusters_after_parallel();
	void initialize_coplanar_clusters_after_parallel();
	void calculate_energy_changing_for_regularization_initialization();
	int get_changed_degrees_of_freedom_after_a_regular_operation(int i);
	void calculate_energy_changing_for_regularization_initialization_normal();
    double euclidean_distance_change_after_a_coplanar_or_parallel(vector<int> this_parallel_cluster, int id_cluster);
    double normal_deviation_change_after_a_coplanar_or_parallel(vector<int> this_parallel_cluster, int id_cluster);
    double normal_deviation_change_after_orthogonal(vector<pair<int, Inexact_Vector_3>> this_orthogonal_cluster);
	double energy_changed_regularization_normal_deviation(double dis, double changed_freedom);
    bool convert_form_regularization_normal(int i, pair<pair<int, vector<int>>, vector<double>>& one_element);
	
    double euclidean_distance_change_after_orthogonal(vector<pair<int, Inexact_Vector_3>>this_orthogonal_cluster);
	

	
	double energy_changed_regularization_euclidean_distance(double dis, double changed_freedom);
	

	Inexact_Point_3 get_specific_centroids(int i);
	void initial_regular_groups_done();
	void update_regular_relations_after_add_remove(int one_change_plane);
    void update_regular_relations_after_splite(int one_splite, vector<int> respective_planes);
	void coset_color();

public:

	double get_x_min() const;
	double get_x_max() const;
	double get_y_min() const;
	double get_y_max() const;
	double get_z_min() const;
	double get_z_max() const;
	double get_bbox_diagonal() const;
	double get_average_spacing() const;
	int get_min_points() const;
	double get_epsilon() const;
	int get_knn() const;
	double get_normal_threshold() const;
	int get_neighborhood_algorithm() const;
	double get_neighborhood_distance() const;
	double get_regularization_lambda() const;
	double get_tolerance_angle() const;
	double get_tolerance_coplanarity() const;

protected:
    string path_point_cloud;
    string path_point_cloud_extension;
    string path_clusters;
	Pwn_vector points;

protected:
    vector<vector<int> > spherical_neighborhood;
    vector<Inexact_Point_3> planes_centroids;

protected:
	double x_min, x_max, y_min, y_max, z_min, z_max;
	double bbox_diagonal;
	double average_spacing;
	bool spacing_is_known;
	int min_points;
	double epsilon;
	int knn;
	bool should_compute_knn;
	double normal_threshold;

	double discretization_angle;
	double discretization_distance;

	int neighborhood_algorithm;
	double neighborhood_distance;
	bool should_compute_neighborhood;
    vector<std::set<int> > neighbors;

	double lambda;
	double tolerance_angle;
	double tolerance_coplanarity;
	

public:

	bool if_oriented_normal;
	int orthogonal_numbers;
	int parallel_numbers;
	int coplanar_numbers;


	int primitives_number;
	double coverage;
	int ori_primitives_number;
	double ori_coverage;
	double ori_mean_error;
	int ori_inliers_number;
	int last_primitives_number;
	double last_normal_deviation;
	double last_coverage;
	double last_mean_error;
	double last_freedom_planes;
    vector<Inexact_Point_3> projected_points;
    vector<Inexact_Plane> planes_0; // original planes
    vector<Inexact_Plane> planes_1; // regularized planes
    vector<Inexact_Plane> planes_2; // discretized planes
	
    vector<vector<double>> planes; // Raphael added this just for debugging
    vector<CGAL::Color> planes_to_colors;
    vector<vector<int> > planes_to_inliers;
    vector<int> inliers_to_planes;
    vector<CGAL::Color> inliers_to_natural_colors;

	int non_coplanar_planes;
    vector<int> coplanarity_table;

    vector<vector<Inexact_Point_3> > alpha_shapes_pts;
    vector<CGAL::Color> alpha_shapes_colors;

    vector<vector<Inexact_Point_3> > best_rectangles_pts;
    vector<CGAL::Color> best_rectangles_colors;

    vector<vector<Inexact_Point_3> > convex_hulls_pts;
    vector<CGAL::Color> convex_hulls_colors;
    vector<pair<int, int> > convex_hulls_seeds;


	bool if_stop;

	bool if_constraint;


	double lambda_r;
	double lambda_regular;
	int weight_mode;
	double lambda_c;
	double lambda_fidelity;

    int metric_type;
    int max_iter;

	double mean_error;
	double all_error;
	double ori_all_error;

	double mean_normal_diviation;
	double mean_distance_diaviation;

	double all_distance_diaviation;
	double all_normal_diaviation;

	size_t number_of_assigned_points;
	double size_current_primitives;
	double mean_distance_current;
	double mean_normal_current;
	double freedom_of_planes_current;

	double ori_mean_normal_diviation;
	double number_inlier_before_opers;
	double ori_freedom_of_planes;
	double freedom_of_planes;


	int old_size_current_primitives;

	double old_mean_distance_diaviation;

	double old_mean_normal_diaviation;
	double old_freedom_of_planes;
	double old_coverage;

    vector<int> if_insert_candidate;
	double mean_distance_diaviation_before_transfer;
	double mean_normal_diviation_befor_transfer;
	int t_m;
	int t_merge;
	int t_split;
	int t_insert;
	int t_exlude;
	int t_regularization;
	int all_t_transfer;
	int all_t_merge;
	int all_t_split;
	int all_t_insert;
	int all_t_exlude;
	int all_regularization;

    vector<int> region_type;

    vector<int> points_if_added;
	int number_of_insert_exclude;
    vector<vector<int>> primitive_connection;

	int number_iterations;
    vector<pair<pair<int, vector<int>>, vector<double>>> good_points_shape;
    vector<pair<pair<int, vector<int>>, vector<double>>> bad_points_shape;
	int t_l;
	double interval_all;
    vector<pair<double, double>> planes_to_unit_sphere_degree;//gauss sphere degree for each plane
    vector<double> planes_to_inliernumber;//inlier number ratio of each plane
    vector<int> planes_to_parallel_clusters;//plane's parallel cluster index

    vector <vector<int>> parallel_clusters_to_planes;//index of planes for each parallel cluester, exist single cluster






    vector<Inexact_Vector_3> normal_parallel_clusters;//normals for parallel clusters
	
	






	
    vector<vector<pair<int, Inexact_Vector_3>>> parallel_id_changed_normal_after_orthogonal;//id of parallel cluster and regularized normal
    vector<Inexact_Point_3> planes_centroids_coplanar;
    vector<vector<vector<int>>> parallel_cluster_to_coplanr_cases;
    vector<vector<double>> energy_changing_initial_regularization;
    vector<bool> if_regularization_can_conducted;
    vector<bool> planes_if_regularized;

    vector<int> planes_to_parallel_done;
    vector<int> planes_to_orthogonal_done;
    vector<int> planes_to_coplanar_done;
    vector<vector<int>> parallel_done_to_planes;
    vector<vector<int>> orthogonal_done_to_planes;
    vector<vector<int>> coplanar_done_to_planes;




};

