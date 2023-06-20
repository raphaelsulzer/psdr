#include "shape_detector_index_map.h"


Shape_Detector_Index_Map::Shape_Detector_Index_Map()
{
}


Shape_Detector_Index_Map::Shape_Detector_Index_Map(std::vector<int> & inliers_to_planes)
	: m_indices (new std::vector<int>(inliers_to_planes.size(), -1))
{
	for (size_t i = 0 ; i < inliers_to_planes.size() ; ++i) {
		(*m_indices)[i] = inliers_to_planes[i];
	}
}


Shape_Detector_Index_Map::~Shape_Detector_Index_Map()
{
}
