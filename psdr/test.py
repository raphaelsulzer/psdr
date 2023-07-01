from pypsdr import psdr
import unittest
import numpy as np
import os

class MyTestCase(unittest.TestCase):

    # def setUp(self):
    #     self.ps = psdr(verbosity=1)
    #
    # def test_a_load_points(self):
    #     self.assertFalse(self.ps.load_points("../example/data/anchor/pointcloud/file.ply"))
    #
    # def test_b_detect(self):
    #     self.assertFalse(self.ps.detect())
    #
    # def test_c_refine(self):
    #     self.assertFalse(self.ps.refine(max_iter=-1))
    #
    # def test_d_save_ply(self):
    #     os.makedirs("../example/data/anchor/convex_hulls",exist_ok=True)
    #     self.assertFalse(self.ps.save("../example/data/anchor/convex_hulls/file.ply",primitive_type="convex"))
    #
    # def test_e_save_npz(self):
    #     os.makedirs("../example/data/anchor/vertex_groups",exist_ok=True)
    #     self.assertFalse(self.ps.save("../example/data/anchor/vertex_groups/file.npz"))

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