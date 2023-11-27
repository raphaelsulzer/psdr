from pypsdr import psdr
import unittest
import numpy as np
import os
# import open3d as o3d

class MyTestCase(unittest.TestCase):

    def test_pipeline(self):
        EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "example","data","anchor")
        print(EXAMPLE_DIR)
        ps = psdr(verbosity=1)
        ps.load_points(os.path.join(os.path.dirname(__file__), "..", "example","anchor.ply"))
        bb_diagonal = ps.get_bounding_box_diagonal()
        # pcd = o3d.io.read_point_cloud(os.path.join(os.path.dirname(__file__), "..", "example","anchor.ply"))
        # bbd = np.linalg.norm(pcd.get_min_bounds(),pcd.get_max_bounds())
        # self.assertEqual(bb_diagonal,bbd)
        ps.detect(epsilon=0.008*bb_diagonal,min_inliers=80)
        ps.refine(max_iterations=-1)
        os.makedirs(os.path.join(EXAMPLE_DIR,"convex_hulls"), exist_ok=True)
        ps.save(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.ply")), primitive_type="convex")
        os.makedirs(os.path.join(EXAMPLE_DIR,"vertex_groups"), exist_ok=True)
        ps.save(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.npz")))
        data = np.load(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.npz")))
        self.assertEqual(data["group_parameters"].shape, (66,4))

if __name__ == '__main__':
    unittest.main()