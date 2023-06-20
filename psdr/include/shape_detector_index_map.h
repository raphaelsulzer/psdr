#pragma once
//#include "defs.h"
#include "defs_cgal_ui.h"
#include <CGAL/Classification/property_maps.h>


class Shape_Detector_Index_Map
{
	boost::shared_ptr<std::vector<int> > m_indices;

public:
	typedef std::size_t key_type;
	typedef int value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	Shape_Detector_Index_Map();

	Shape_Detector_Index_Map(std::vector<int> & inliers_to_planes);

	~Shape_Detector_Index_Map();

	inline friend value_type get(const Shape_Detector_Index_Map & M, const key_type & K)
	{
		return (*(M.m_indices))[K];
	}
};

