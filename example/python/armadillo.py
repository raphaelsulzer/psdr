import os
from pypsdr import psdr


params = dict()
params["40"] = {"min_inliers":4000,"epsilon":5,"normal_th":0.5}
params["100"] = {"min_inliers":1000,"epsilon":4,"normal_th":0.5}
params["500"] = {"min_inliers":120,"epsilon":0.9,"normal_th":0.75}
params["2000"] = {"min_inliers":50,"epsilon":0.25,"normal_th":0.75}
params["10000"] = {"min_inliers":25,"epsilon":0.06,"normal_th":0.75}

model = "armadillo"

for i,p in params.items():

    ps = psdr(verbosity=1)

    file = "../data/{}/pointcloud/file.ply".format(model)
    os.makedirs(os.path.dirname(file),exist_ok=True)

    ps.load_points(file)

    ps.detect(**p)
    ps.refine(max_iter=20)

    file = "../data/{}/refined/file{}.npz".format(model,i)
    os.makedirs(os.path.dirname(file),exist_ok=True)
    ps.save(file)

    file = "../data/{}/alpha/file{}.ply".format(model,i)
    os.makedirs(os.path.dirname(file),exist_ok=True)
    ps.save(file,"alpha")

    file = "../data/{}/pointgroups/file{}.ply".format(model,i)
    os.makedirs(os.path.dirname(file),exist_ok=True)
    ps.save(file,"pointcloud")


