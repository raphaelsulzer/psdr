# PSDR: Planar shape detection and refinement

![build](https://github.com/raphaelsulzer/psdr/actions/workflows/build.yml/badge.svg)

This repository contains a pipeline for planar shape detection [[1]](#references) and refinement [[2]](#references) from point clouds. 
The source code is written in C++ and Python bindings are provided for the main functionality.

<p float="left">
  <img style="width:800px;" src="./media/city.gif">
</p>

# TODO: make a build.yml and test.yml for in github actions which builds and tests the repo at every push


## :clipboard: Features

- Reading of point clouds (.ply) or vertex groups ([.vg](https://abspy.readthedocs.io/en/latest/vertexgroup.html), .npz) as input
- Planar shape detection based on a robust and efficient region growing algorithm [[1]](#references) (also in [CGAL](https://doc.cgal.org/latest/Shape_detection/index.html#Shape_detection_RegionGrowing))
- Planar shape refinement based on an optimization that seeks the best trade-off between fidelity, completeness and simplicity of the configuration [[2]](#references)
- Writing of planar shapes as 2D convex hulls, alpha shapes or minimal rectangles (.ply) or as vertex groups ([.vg](https://abspy.readthedocs.io/en/latest/vertexgroup.html), .npz ,.ply).

## :bricks: Installation

First, clone this repository, create a new conda environment called psdr and install all necessary dependencies. 
Then install PSDR including its Python bindings  with `pip install .`. Make sure that the conda environment is activated for all steps.
If all tests complete successfully you are ready to use PSDR.

```
git clone https://github.com/raphaelsulzer/psdr.git
cd psdr/psdr
conda create --name psdr
conda activate psdr
conda install -y -c conda-forge xtensor xtensor-io spdlog cgal anaconda::mpfr yaml-cpp omnia::eigen3
pip install .
python -m unittest test.py               
```



## :computer: Usage


### Python

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
ps.save(example/data/anchor/point_groups.ply,"pointcloud")               
ps.save(example/data/anchor/point_groups.vg)                              
ps.save(example/data/anchor/point_groups.npz)                              
```
For more Python examples see `example/python`.

### C++


```
auto SD = Shape_Detector();
SD.load_points(example/data/anchor/pointcloud.ply);
SD.set_detection_parameters(20,0.02,0.8,10);
auto SC = Shape_Container(&SD);
SC.detect();
SC.refine(10);
SC.save(example/data/gargoyle/groups.npz);
SC.save(example/data/gargoyle/rectangles.ply,"rectangle");
```
For a cmake project that uses PSDR see `example/cpp`.



## :camera_flash: Examples

### Refinement
<p float="left">
  <img style="width:800px;" src="./media/anchor.png">
</p>

### Levels of detail
<p float="left">
  <img style="width:800px;" src="./media/armadillo.png">
</p>

### Representations
<p float="left">
  <img style="width:800px;" src="./media/tarbosaurus.png">
</p>




## :book: References

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

## Citation

If you use this library in your work, please consider citing the above papers and this repository:
```bibtex
@misc{sulzer2023psdr,
  author = {Sulzer, Raphael},
  title = {PSDR: Planar shape detection and refinement},
  year = {2023},
  howpublished = {GitHub Repository},
  url = {https://github.com/raphaelsulzer/psdr}
}
```

