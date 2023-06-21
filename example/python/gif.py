import numpy as np
from PIL import Image
from glob import glob
import os
from copy import deepcopy

images = []


path = "/home/rsulzer/cpp/psdr/example/data/"
model = "gargoyle"


files = glob(os.path.join(path,model,'renders','*'))
files.sort(key=os.path.getmtime)

for f in files:

    image = Image.open(os.path.join(path,model,'renders',f))

    new_image = Image.new("RGBA", image.size, "WHITE")  # Create a white rgba background
    new_image.paste(image, (0, 0), image)  # Paste the image on the background. Go to the links given below for details.
    new_image.convert('RGB').save('test.jpg', "JPEG")

    images.append(new_image)




images = np.array(images)
images = np.hstack((images,np.flip(images)))

images = list(images)

images[0].save(os.path.join(path, model + ".gif"),
               save_all=True, append_images=images[1:], optimize=True, duration=100, loop=0)


