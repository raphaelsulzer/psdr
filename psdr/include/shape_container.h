#pragma once
#include "defs_cgal_ui.h"
#include "shape_detector.h"


class Shape_Container
{
public:
	Shape_Container();

	Shape_Container(Shape_Detector* _SD);

    int detect();

	void regularize();

    int refine(int max_iter = -1, int max_seconds = -1);
	void copy_primitives_from_detector();
	
	void copy_support_planes_from_detector();

	void discard_degenerate_primitives();

	bool is_empty();

	int get_number_of_primitives();

	int get_number_of_support_planes();

	size_t size_of_inexact_convex_hulls() const;

	size_t size_of_inexact_best_rectangles() const;

	const std::vector<Inexact_Point_3> & get_triangle(const int triangle_index) const;

	std::pair<int, int> get_range_of_triangles(const size_t primitive_index) const;

	const std::vector<Inexact_Point_3> & get_convex_hull_points(const size_t primitive_index) const;

	CGAL::Color get_primitive_color(const size_t primitive_index) const;

	const std::vector<Inexact_Point_3> & get_best_rectangle_points(const size_t primitive_index) const;

	bool is_primitive_degenerate(const size_t primitive_index) const;

    void get_inliers(size_t primitive_index, vector<Point_with_normal> & inliers);
    void get_inliers(size_t primitive_index, vector<size_t> & inliers);

	size_t get_total_number_of_inliers();

	const Inexact_Plane & get_plane(const size_t p) const;

public:

	Shape_Detector* SD;

    void to_ply(const string& filename, const string& type="convex");
    void inliers_to_ply(const string & filename);
    int save(const string & filename, const string& type="convex");


protected:
	std::vector<std::vector<Inexact_Point_3> > inexact_alpha_shapes;
	std::vector<std::vector<Inexact_Point_3> > inexact_convex_hulls;
	std::vector<std::pair<int, int> > convex_hulls_to_alpha_shapes;
	std::vector<std::vector<Inexact_Point_3> > inexact_best_rectangles;

	int degenerate_primitives;
	std::vector<bool> is_not_degenerate;
	
	std::vector<Inexact_Plane> support_planes_of_primitives;
};

