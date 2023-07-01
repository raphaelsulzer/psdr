import os
from pypsdr import psdr


params = dict()
params["anchor"] = {"epsilon":0.10,"min_inliers":50,"knn":10,"normal_th":0.8}

model = "anchor"


ps = psdr(verbosity=1)

file = "../data/{}/pointcloud/file.ply".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)

# file = "../data/{}/convexes_detected/file.vg".format(model)

ps.load_points(file)

ps.detect(**params[model])

file = "../data/{}/convexes_detected/file.ply".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)
ps.save(file)
file = "../data/{}/pointgroups_detected/file.ply".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)
ps.save(file,"pointcloud")

ps.refine()

file = "../data/{}/convexes_refined/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file)
file = "../data/{}/pointgroups_refined/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file,"pointcloud")