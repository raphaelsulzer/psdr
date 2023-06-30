import numpy as np
from PIL import Image
from glob import glob
import os
from copy import deepcopy
import subprocess
import shutil
import numpy as np
from tqdm import tqdm

path = "/home/rsulzer/cpp/psdr/example/data/"
# path = "/home/rsulzer/python/compod/example/data"
model = "city"
mode1 = "pointcloud"
mode2 = "convexes_refined"
# mode2 = None

files = glob(os.path.join(path,model,'renders',mode1,'*'))
files.sort(key=os.path.getmtime)

os.makedirs(os.path.join(path,model,'renders',mode1+"_"),exist_ok=True)


def get_concat_h(im1, im2):
    dst = Image.new('RGBA', (im1.width + im2.width, im1.height),"WHITE")
    dst.paste(im1, (0, 0), im1)
    dst.paste(im2, (im1.width, 0), im2)
    return dst

def get_concat_v(im1, im2):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst



if True:

    pause = 100
    show_edges = np.ones(pause,dtype=int)
    show_edges[np.random.randint(0,pause-1,30)] = 0

    # os.makedirs(os.path.join(path,model,'rendersw'),exist_ok=True)
    for i,f in enumerate(tqdm(files)):

        f = os.path.split(f)[-1]

        image1 = Image.open(os.path.join(path,model,'renders',mode1,f))

        if mode2 is not None:
            image2 = Image.open(os.path.join(path,model,'renders',mode2,f))
            new_image = get_concat_h(image1,image2)
        else:
            new_image = image1

        # new_image = Image.new("RGBA", image1.size, "WHITE")  # Create a white rgba background
        # new_image.paste(image1, (0, 0), image1)

        i1 = i
        i2 = 2*len(files)-i-1+pause

        new_image = new_image.convert('RGB')

        new_image.save(os.path.join(path,model,'renders',mode1+"_",str(i1)+".png"), "png")
        new_image.save(os.path.join(path,model,'renders',mode1+"_",str(i2)+".png"), "png")

    for i in range(pause):

        dst = os.path.join(os.path.join(path,model,'renders',mode1+"_",str(i1+i+1)+".png"))
        if show_edges[i]:
            src = os.path.join(path,model,'renders',mode1+"_",str(i1)+".png")
            shutil.copyfile(src, dst)
        else:
            src = os.path.join(path,model,'renders',mode1+"_",str(i1-1)+".png")
            shutil.copyfile(src, dst)




# TODO: save outputs to media folder, first frame as a still and the gif

outfile = str(os.path.join(path,'..','..','media',"{}.gif".format(mode1)))
os.makedirs(os.path.split(outfile)[0],exist_ok=True)

command = '/root/Downloads/gifski-1.11.0/linux/gifski -r 20 --quality 100 --motion-quality 100 --lossy-quality 100 --output '\
+outfile+" "+str(os.path.join(path, model, "renders",mode1+"_",'*.png'))
command += " --width 800"
print(*command)
p=subprocess.Popen(command,shell=True)
p.wait()
