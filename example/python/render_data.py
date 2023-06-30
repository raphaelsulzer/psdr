import math

model_dict = dict()
model_dict["city"] = dict()
model_dict["city"]["rotation"] = None
model_dict["city"]["light"] = 20000
model_dict["city"]["point_size"] = 0.03
model_dict["city"]["light_rotation"] = [math.pi / 5, -math.pi / 5, 0]

model_dict["bunny"] = dict()
model_dict["bunny"]["rotation"] = [-math.pi / 2, 0, 0]
model_dict["bunny"]["light"] = 3 ** 3
model_dict["bunny"]["point_size"] = 0.001
model_dict["bunny"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]

model_dict["anchor"] = dict()
model_dict["anchor"]["rotation"] = None
model_dict["anchor"]["light"] = 50000
model_dict["anchor"]["point_size"] = 0.4
model_dict["anchor"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]

model_dict["armadillo"] = dict()
model_dict["armadillo"]["rotation"] = None
model_dict["armadillo"]["light"] = 50000
model_dict["armadillo"]["point_size"] = 0.4
model_dict["armadillo"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]