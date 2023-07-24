import os
from pypsdr import psdr


params = dict()
params["city"] = {"min_inliers":120,"epsilon":0.05,"knn":10,"normal_th":0.9}

model = "city"

ps = psdr(verbosity=1)

file = "../data/{}/pointcloud/file.ply".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)


ps.load_points(file)

ps.detect(**params[model])

file = "../data/{}/convexes_detected/file.ply".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)
ps.save(file)
file = "../data/{}/convexes_detected/file.npz".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)
ps.save(file)

ps.refine()

file = "../data/{}/convexes_refined/file.ply".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file)
file = "../data/{}/convexes_refined/file.npz".format(model)
os.makedirs(os.path.dirname(file), exist_ok=True)
ps.save(file)