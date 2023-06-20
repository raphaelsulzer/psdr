# PSDR: Planar shape detection and refinement

TODO: rotating gif 

This repository contains a pipeline for planar shape detection [1] and refinement [2] from point clouds. The source code is written in C++ and we also provide convenient Python bindings.

# Features


- Reading of point clouds (.ply) or vertex groups ([.vg](https://abspy.readthedocs.io/en/latest/vertexgroup.html), .npz) as input
- Planar shape detection based on a robust and efficient region growing algorithm [1] (also implemented in CGAL)
- Planar shape refinement based on an optimization finding the best trade-off between fidelity, completeness and simplicity of the configuration [2]
- Writing of planar shapes as convex hulls, alpha shapes or minimal rectangles (.ply) or as vertex groups ([.vg](https://abspy.readthedocs.io/en/latest/vertexgroup.html), .npz).

# Installation

Simply clone the repository and install in a new conda environment using pip:

```
git clone https://github.com/raphaelsulzer/psdr.git
cd psdr
conda create --name psdr
conda activate psdr
pip install . 
```

You are now ready to use PSDR.


# Usage


## Python

```
from pypsdr import psdr

# initialise a planar shape detector and load input points                                              
ps = psdr(verbosity=1)                                               
ps.load_points(example/data/anchor/pointcloud.ply)

# detect planar shapes with default values
ps.detect(min_inliers=20,epsilon=0.02,normal_th=0.8,knn=10)

# refine planar shape configuration until convergence (i.e. no limit on number of iterations)
ps.refine(max_iter=-1)

# export planar shapes and vertex groups  
ps.save(example/data/anchor/convexes.ply,"convex")                  
ps.save(example/data/anchor/rectangles.ply,"rectangles")            
ps.save(example/data/anchor/alpha_shapes.ply,"alpha")               
ps.save(example/data/anchor/groups.vg)                              
ps.save(example/data/anchor/groups.npz)                             
```


## C++

See `example/cpp` for a full example of a project that uses PSDR.

```
auto SD = Shape_Detector();
SD.load_points();
SD.set_detection_parameters(min_inliers, pd_epsilon, knn, normal_th);
auto SC = Shape_Container(&SD);
SC.detect(20,0.02,0.8,10);
SC.refine();
SC.save(example/data/anchor/groups.npz);
```



# Examples

City 
point cloud - convexes with random color - refined convex shapes

Angel 
convexes with decreasing epsilon and num_inliers

Gargoyle with rectangle, alpha shape and convex hulls


# References

```bibtex
@article{1,
  title={Creating large-scale city models from 3D-point clouds: a robust approach with hybrid representation},
  author={Lafarge, Florent and Mallet, Cl{\'e}ment},
  journal={International journal of computer vision},
  volume={99},
  pages={69--85},
  year={2012},
  publisher={Springer}
}
```

```bibtex
@inproceedings{2,
  title={Finding Good Configurations of Planar Primitives in Unorganized Point Clouds},
  author={Yu, Mulin and Lafarge, Florent},
  booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},
  pages={6367--6376},
  year={2022}
}
```