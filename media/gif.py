import numpy as np
from PIL import Image, ImageDraw, ImageFont
from glob import glob
import os
from copy import deepcopy
import subprocess
import shutil
import numpy as np
from tqdm import tqdm


def get_concat_h(images):
    im1 = images[0]
    im2 = images[1]

    dst = Image.new('RGBA', (im1.width + im2.width, im1.height),"WHITE")
    dst.paste(im1, (0, 0), im1)
    dst.paste(im2, (im1.width, 0), im2)
    return dst

def get_concat_v(images):
    im1 = images[0]
    im2 = images[1]

    dst = Image.new('RGBA', (im1.width, im1.height + im2.height),"WHITE")
    dst.paste(im1, (0, 0), im1)
    dst.paste(im2, (0, im1.height), im2)
    return dst

def get_concat_4(images):
    im1 = images[0]
    im2 = images[1]
    im3 = images[2]
    im4 = images[3]

    dst = Image.new('RGBA', (2*im1.width, 2*im1.height),"WHITE")
    dst.paste(im1, (0, 0), im1)
    dst.paste(im2, (im1.width, 0), im2)
    dst.paste(im3, (0, im1.height), im3)
    dst.paste(im4, (im1.width, im1.height), im4)
    return dst

def draw_text(im,text):

    font = ImageFont.truetype("/usr/share/fonts/timr45w.ttf", 24)

    if '\n' in text:
        text_size = font.getsize_multiline(text)
    else:
        text_size = font.getsize(text)

    margin = 20
    # set button size + 10px margins
    button_size = (text_size[0] + margin, text_size[1] + margin)

    # create image with correct size and black background
    button_img = Image.new('RGBA', button_size, (255,255,255,200))

    # put text on button with 10px margins
    button_draw = ImageDraw.Draw(button_img)
    button_draw.multiline_text((margin/2, margin/2), text, (0,0,0), font=font)

    im.alpha_composite(button_img,(margin,im.width-margin-button_size[1]))
    # im.paste(button_img,(margin,im.width-margin-button_size[1]),button_img)

    # draw = ImageDraw.Draw(im)
    # # font = ImageFont.truetype(<font-file>, <font-size>)
    # # to make any font work search for the otf or ttf file: https://pillow.readthedocs.io/en/stable/reference/ImageFont.html#PIL.ImageFont.truetype
    # # font = ImageFont.load_default(25)
    # # draw.text((x, y),"Sample Text",(r,g,b))
    # draw.text((im.width / 2, im.height / 2), text, (0, 0, 0), font=font)

def prep_images(modes,out_path,text=None,slow_start=50,pause=100, arrange=None):

    default_path = "/home/rsulzer/cpp/psdr/example/data/"

    if not isinstance(modes,list):
        modes = [modes]

    files = glob(os.path.join(path_dict.get(modes[0],default_path), model, 'renders', modes[0], '*'))
    files.sort(key=os.path.getmtime)

    for ss in range(slow_start):
        files = [files[0]]+files

    # files = files[:10]

    os.makedirs(os.path.join(path_dict.get(modes[0],default_path), model, 'renders', modes[0] + "_"), exist_ok=True)

    if pause > 0:
        show_edges = np.ones(pause,dtype=int)
        show_edges[np.random.randint(0,pause-1,30)] = 0
        show_edges[int(pause/2):] = 1

    # os.makedirs(os.path.join(path,model,'rendersw'),exist_ok=True)
    for i,f in enumerate(tqdm(files)):

        f = os.path.split(f)[-1]

        images = []

        for j,mode in enumerate(modes):
            im = Image.open(os.path.join(path_dict.get(modes[0],default_path), model, 'renders', mode, f))
            images.append(im)

            if text is not None:
                draw_text(im,text[j])


        if arrange == "h":
            new_image = get_concat_h(images)
        elif arrange == "v":
            new_image = get_concat_v(images)
        elif arrange == "4":
            new_image = get_concat_4(images)
        else:
            new_image = Image.new("RGBA", images[0].size, "WHITE")  # Create a white rgba background
            new_image.paste(images[0], (0, 0), images[0])

        i1 = i
        i2 = 2*len(files)-i-1+pause

        # new_image = new_image.convert('RGB')

        new_image.save(os.path.join(out_path,str(i1)+".png"), "png")
        new_image.save(os.path.join(out_path,str(i2)+".png"), "png")

    for i in range(pause):

        dst = os.path.join(os.path.join(out_path,str(i1+i+1)+".png"))
        if show_edges[i]:
            src = os.path.join(out_path,str(i1)+".png")
            shutil.copyfile(src, dst)
        else:
            src = os.path.join(out_path,str(i1-1)+".png")
            shutil.copyfile(src, dst)



def make_gif(in_file,out_file):

    os.makedirs(os.path.split(out_file)[0],exist_ok=True)

    # 24 to 30 fps is real time

    command = '/usr/share/gifski-1.11.0/linux/gifski --fps 24 --quality 100 --motion-quality 100 --lossy-quality 100 --output '+out_file+" "+in_file
    command += " --width 800"
    print(*command)
    p=subprocess.Popen(command,shell=True)
    p.wait()


if __name__ == "__main__":

    path_dict = dict()
    path_dict["pointcloud"] = "/home/rsulzer/cpp/psdr/example/data/"
    path_dict["convexes_refined"] = "/home/rsulzer/cpp/psdr/example/data/"
    path_dict["dense_mesh"] = "/home/rsulzer/python/compod/example/data"
    path_dict["colored_soup"] = "/home/rsulzer/python/compod/example/data"
    path_dict["polygon_mesh"] = "/home/rsulzer/python/compod/example/data"





    ## city (compod)
    # model = "city"
    # text = ["Input polygons","Polygon mesh (COMPOD)","Triangle mesh (COMPOD)","Triangle mesh"]
    # methods = ["convexes_refined","colored_soup","polygon_mesh","dense_mesh"]

    ## city (psdr)
    model = "city"
    # methods = ["pointcloud", "convexes_refined"]
    # text = ["Input pointcloud", "Detected shapes"]

    ## armadillo
    # model = "armadillo"
    # methods = ["alpha100","alpha500","alpha2000","alpha10000"]
    # text = ["epsilon = 4\n#inliers = 1000",
    #         "epsilon = 0.9\n#inliers = 120",
    #         "epsilon = 0.25\n#inliers = 50",
    #         "epsilon = 0.06\n#inliers = 25",]

    # ## tarbosaurus
    # model = "tarbosaurus"
    # text = ["Vertex groups","Convex hulls","Alpha shapes","Minimal rectangles"]
    # methods = ["pointgroups","convex_hull","alpha_shape","rectangle"]

    ## anchor
    # model = "anchor"
    # text = ["Detected","Refined","Detected","Refined"]
    # methods = ["pointgroups_detected","pointgroups_refined","convexes_detected","convexes_refined"]
    #

    # outfpath for images and outfile for gif  PSDR
    # out_path = "/home/rsulzer/cpp/psdr/example/data/{}/renders/gif".format(model)
    # out_file = "/home/rsulzer/cpp/psdr/media/{}.gif".format(model)

    # outfpath for images and outfile for gif COMPOD
    out_path = "/home/rsulzer/python/compod/example/data/{}/renders/gif".format(model)
    out_file = "/home/rsulzer/python/compod/media/{}.gif".format(model)


    # os.makedirs(out_path,exist_ok=True)
    #
    # prep_images(methods,out_path,text=text,slow_start=0,pause=0,arrange='4')
    # prep_images(["convexes_refined"],out_path,pause=80)

    # out_file = "home/rsulzer/cpp/psdr/media/{}.gif".format(model)
    os.makedirs(os.path.split(out_file)[0],exist_ok=True)
    make_gif(out_path+"/*.png",out_file)



