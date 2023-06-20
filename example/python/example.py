from pypsdr import psdr


# anchor

ps = psdr(verbosity=1)
ps.load_points("../data/anchor/anchor_1.ply")
ps.detect(epsilon=0.2)
ps.save("../data/anchor/convexes_detected.ply")
ps.refine()
ps.save("../data/anchor/convexes_refined.ply")
ps.save("../data/anchor/alpha-shapes_refined.ply","alpha")


# angel

