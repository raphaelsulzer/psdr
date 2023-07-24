from pypsdr import psdr
import unittest
import numpy as np
import os

class MyTestCase(unittest.TestCase):



    def test_pipeline(self):
        EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "example","data","anchor")
        print(EXAMPLE_DIR)
        ps = psdr(verbosity=1)
        ps.load_points(os.path.join(EXAMPLE_DIR,"pointcloud","file.ply"))
        ps.detect(epsilon=0.10,min_inliers=50)
        ps.refine(max_iter=-1)
        os.makedirs(os.path.join(EXAMPLE_DIR,"convex_hulls"), exist_ok=True)
        ps.save(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.ply")), primitive_type="convex")
        os.makedirs(os.path.join(EXAMPLE_DIR,"vertex_groups"), exist_ok=True)
        ps.save(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.npz")))
        data = np.load(str(os.path.join(EXAMPLE_DIR,"convex_hulls","file.npz")))
        self.assertEqual(data["group_parameters"].shape, (130,4))

if __name__ == '__main__':
    unittest.main()