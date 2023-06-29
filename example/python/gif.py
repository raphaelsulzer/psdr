import numpy as np
from PIL import Image
from glob import glob
import os
from copy import deepcopy
import subprocess

# images = []

path = "/home/rsulzer/cpp/psdr/example/data/"
path = "/home/rsulzer/python/compod/example/data"
model = "city"
mode = "polygon_mesh"

files = glob(os.path.join(path,model,'renders',mode,'*'))
files.sort(key=os.path.getmtime)

os.makedirs(os.path.join(path,model,'renders',mode+"_"),exist_ok=True)

if False:

    # os.makedirs(os.path.join(path,model,'rendersw'),exist_ok=True)
    for i,f in enumerate(files):

        image = Image.open(os.path.join(path,model,'renders',mode,f))

        new_image = Image.new("RGBA", image.size, "WHITE")  # Create a white rgba background
        new_image.paste(image, (0, 0), image)

        i1 = i
        i2 = 2*len(files)-i-1

        new_image = new_image.convert('RGB')

        new_image.save(os.path.join(path,model,'renders',mode+"_",str(i1)+".png"), "png")
        new_image.save(os.path.join(path,model,'renders',mode+"_",str(i2)+".png"), "png")




command = ['ffmpeg','-framerate','10','-y','-i',os.path.join(path, model, "renders",mode+"_","%d.png"), os.path.join(path,model,'renders',"{}.gif".format(mode))]
print(*command)
p=subprocess.Popen(command)
p.wait()
