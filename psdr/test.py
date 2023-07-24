from pypsdr import psdr
import unittest
import numpy as np
import os

class MyTestCase(unittest.TestCase):



    def test_pipeline(self):
        ps = psdr(verbosity=1)
        ps.load_points("../example/data/anchor/pointcloud/file.ply")
        ps.detect()
        ps.refine(max_iter=-1)
        os.makedirs("../example/data/anchor/convex_hulls", exist_ok=True)
        ps.save("../example/data/anchor/convex_hulls/file.ply", primitive_type="convex")
        os.makedirs("../example/data/anchor/vertex_groups", exist_ok=True)
        ps.save("../example/data/anchor/vertex_groups/file.npz")
        data = np.load("../example/data/anchor/vertex_groups/file.npz")
        self.assertEqual(data["group_parameters"].shape, (133,4))

if __name__ == '__main__':
    unittest.main()