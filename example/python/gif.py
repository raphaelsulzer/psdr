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
path = "/home/rsulzer/python/compod/example/data"
model = "city"
mode = "polygon_mesh"

files = glob(os.path.join(path,model,'renders',mode,'*'))
files.sort(key=os.path.getmtime)

os.makedirs(os.path.join(path,model,'renders',mode+"_"),exist_ok=True)



if True:

    pause = 100
    show_edges = np.ones(pause,dtype=int)
    show_edges[np.random.randint(0,pause-1,30)] = 0

    # os.makedirs(os.path.join(path,model,'rendersw'),exist_ok=True)
    for i,f in enumerate(tqdm(files)):

        image = Image.open(os.path.join(path,model,'renders',mode,f))

        new_image = Image.new("RGBA", image.size, "WHITE")  # Create a white rgba background
        new_image.paste(image, (0, 0), image)

        i1 = i
        i2 = 2*len(files)-i-1+pause


        new_image = new_image.convert('RGB')

        new_image.save(os.path.join(path,model,'renders',mode+"_",str(i1)+".png"), "png")
        new_image.save(os.path.join(path,model,'renders',mode+"_",str(i2)+".png"), "png")

    for i in range(pause):

        dst = os.path.join(os.path.join(path,model,'renders',mode+"_",str(i1+i+1)+".png"))
        if show_edges[i]:
            src = os.path.join(path,model,'renders',mode+"_",str(i1)+".png")
            shutil.copyfile(src, dst)
        else:
            src = os.path.join(path,model,'renders',mode+"_",str(i1-1)+".png")
            shutil.copyfile(src, dst)



# command = ['ffmpeg','-framerate','20','-y','-i',os.path.join(path, model, "renders",mode+"_","%d.png"), os.path.join(path,model,'renders',"{}.gif".format(mode))]
# command = ['/root/Downloads/gifski-1.11.0/linux/gifski','-r','20',
#            # '--quality','100','--motion-quality','100','--lossy-quality','100',
#            '--output',os.path.join(path,model,'renders',"{}.gif".format(mode)),
#            os.path.join(path, model, "renders",mode+"_",'*.png')]
command = '/root/Downloads/gifski-1.11.0/linux/gifski -r 20 --quality 100 --motion-quality 100 --lossy-quality 100 --output '\
+str(os.path.join(path,model,'renders',"{}.gif".format(mode)))+" "+str(os.path.join(path, model, "renders",mode+"_",'*.png'))
print(*command)
p=subprocess.Popen(command,shell=True)
p.wait()
