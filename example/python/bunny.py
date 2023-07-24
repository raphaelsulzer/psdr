import os
from pypsdr import psdr
from compod import all

ps = psdr(verbosity=1)

model = "bunny"

file = "../data/{}/convexes_detected/file.vg".format(model)
os.makedirs(os.path.dirname(file),exist_ok=True)

ps.load_points(file)
ps.detect()

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




