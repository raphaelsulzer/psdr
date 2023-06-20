#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/Color.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h> 
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h> 
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Fuzzy_sphere.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>


typedef CGAL::Lazy_exact_nt<CGAL::Epeck_ft> ET;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_2 Exact_Point_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_3 Exact_Point_3;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Segment_2 Exact_Segment_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Segment_2 Exact_Segment_3;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Line_2 Exact_Line_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Line_3 Exact_Line_3;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Vector_2 Exact_Vector_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Vector_3 Exact_Vector_3;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Plane_3 Exact_Plane;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Direction_2 Exact_Direction_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Triangle_2 Exact_Triangle_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel::Triangle_3 Exact_Triangle_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 Inexact_Point_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Inexact_Point_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Vector_2 Inexact_Vector_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Vector_3 Inexact_Vector_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 Inexact_Plane;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Triangle_3 Inexact_Triangle_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Triangle_2 Inexact_Triangle_2;

// CGAL::Tag_true needs to be passed according to this https://github.com/CGAL/cgal/issues/1700
// without it, assertions are thrown
typedef CGAL::Alpha_shape_vertex_base_2<CGAL::Exact_predicates_inexact_constructions_kernel,CGAL::Default,CGAL::Tag_true>    Vb2;
typedef CGAL::Alpha_shape_face_base_2<CGAL::Exact_predicates_inexact_constructions_kernel,CGAL::Default,CGAL::Tag_true>    Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2> Tds2; 
typedef CGAL::Delaunay_triangulation_2<CGAL::Exact_predicates_inexact_constructions_kernel, Tds2> Triangulation_2; 
typedef CGAL::Alpha_shape_2<Triangulation_2,CGAL::Tag_true>  Alpha_Shape;
typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL::Exact_predicates_inexact_constructions_kernel> CDT2;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT2> Criteria; 
typedef CGAL::Delaunay_mesher_2<CDT2, Criteria> Meshing_engine;

typedef CGAL::Search_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> SearchTraits_3;
typedef CGAL::Kd_tree<SearchTraits_3> Tree_radius;
typedef CGAL::Fuzzy_sphere<SearchTraits_3> Fuzzy_sphere;
typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits_3> Neighbor_search;
typedef Neighbor_search::Tree Tree;

typedef CGAL::Color Color;

typedef CGAL::Surface_mesh<Inexact_Point_3>   Inexact_Mesh;

typedef  Inexact_Mesh::Face_range      Face_range;
typedef  Face_range::iterator       Face_range_iter;
typedef  Inexact_Mesh::Halfedge_index      Halfedge_descriptor;

typedef CGAL::Cartesian_converter<CGAL::Exact_predicates_inexact_constructions_kernel, CGAL::Exact_predicates_exact_constructions_kernel> Ix_2_ex;
typedef CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, CGAL::Exact_predicates_inexact_constructions_kernel> Ex_2_ix;

struct FaceInfo2
{
	FaceInfo2() {}
	int nesting_level;
	bool in_domain() {
		return nesting_level % 2 == 1;
	}
};
typedef CGAL::Triangulation_vertex_base_2<CGAL::Exact_predicates_inexact_constructions_kernel>                      Vb_simplify;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, CGAL::Exact_predicates_inexact_constructions_kernel>    Fbb_simplify;
typedef CGAL::Constrained_triangulation_face_base_2<CGAL::Exact_predicates_inexact_constructions_kernel, Fbb_simplify>        Fb_simplify;
typedef CGAL::Triangulation_data_structure_2<Vb_simplify, Fb_simplify>               TDS_simplify;
typedef CGAL::Exact_predicates_tag                                Itag_simplify;
typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL::Exact_predicates_inexact_constructions_kernel, TDS_simplify, Itag_simplify>  CDT_simplify;
typedef CDT_simplify::Point                                                Point_simplify;
typedef CGAL::Polygon_2<CGAL::Exact_predicates_inexact_constructions_kernel>                                        Polygon_2_simplify;
typedef CDT_simplify::Face_handle                                          Face_handle_simplify;

typedef boost::graph_traits<Inexact_Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Inexact_Mesh>::face_descriptor   face_descriptor;
typedef boost::graph_traits<Inexact_Mesh>::edge_descriptor   edge_descriptor;
