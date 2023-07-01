import os
from pypsdr import psdr

params = dict()
params["tarbosaurus"] = {"min_inliers":100,"epsilon":0.6,"knn":10,"normal_th":0.8}

model = "tarbosaurus"

ps = psdr(verbosity=1)

file = "../data/{}/pointcloud/file.ply".format(model)

ps.load_points(file)
ps.detect(**params[model])
ps.refine()

file = "../data/{}/pointgroups/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file,"pointcloud")

file = "../data/{}/convex_hull/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file,"convex")

file = "../data/{}/alpha_shape/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file,"alpha")

file = "../data/{}/rectangle/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file,"rectangle")