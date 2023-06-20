import math
import pydevd_pycharm
# pydevd_pycharm.settrace('localhost', port=1090, stdoutToServer=True, stderrToServer=True)
import os

import bpy
import numpy as np
from pathlib import Path

import blender_plots as bplt
from scipy.spatial.transform import Rotation

from glob import glob
from pyntcloud import PyntCloud

import sys

from mathutils import Vector, Matrix

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm


imtype = "png"

class MplColorHelper:
    """
    great code from here: https://stackoverflow.com/a/26109298
    which given two values start_val, stop_val makes a color gradient cmap in between.
    then an array passed to cmap.get_rgb gives the corresponding colors in rgb
    values of the array outside [start_val, stop_val] are simply assigned the endoint of the color gradient
    works great in combination with start_val = np.percentile(array,5) and stop_val = np.percentile(array,95)
    """

    def __init__(self, cmap_name, start_val, stop_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
        return self.scalarMap.to_rgba(val)


class BlenderRender:

    def __init__(self):

        self.remove_model = True

        self.rotate = False

        collection = bpy.data.collections.get("Collection")
        if collection:
            for obj in collection.objects:
                bpy.data.objects.remove(obj, do_unlink=True)

            bpy.data.collections.remove(collection)

        ## if collection already exists, remove it and all its objects
        collection = bpy.data.collections.get("MyCollection")
        if collection:
            for obj in collection.objects:
                bpy.data.objects.remove(obj, do_unlink=True)

            bpy.data.collections.remove(collection)

        ## New Collection
        self.coll = bpy.data.collections.new("MyCollection")

        ## Add collection to scene collection
        bpy.context.scene.collection.children.link(self.coll)

        self.scene_coll = bpy.context.scene.collection




    def add_camera(self, cam_file, resolution=(1024, 1024)):



        cam_data = np.load(cam_file)
        # get with C.scene.camera.location
        location = cam_data["location"]
        # get with C.scene.camera.matrix_world.to_euler()
        orientation = cam_data["orientation"]

        self.light_location = location


        ## make camera and link it
        camera_data = bpy.data.cameras.new("Camera")
        self.camera = bpy.data.objects.new("Camera", camera_data)

        # get camera location with C.scene.camera.location
        self.camera.location = location
        # get camera angle with C.scene.camera.matrix_world.to_euler()
        self.camera.rotation_euler = orientation
        self.coll.objects.link(self.camera)


        # change camera size
        bpy.context.scene.render.resolution_x = resolution[0]
        bpy.context.scene.render.resolution_y = resolution[1]

        bpy.context.scene.camera = self.camera



    def add_light(self, location, energy=15):

        # create light datablock, set attributes
        light_data = bpy.data.lights.new(name="light", type='POINT')
        light_data.energy = energy

        # create new object with our light datablock
        self.light = bpy.data.objects.new(name="light", object_data=light_data)

        # link light object
        self.coll.objects.link(self.light)

        # make it active
        bpy.context.view_layer.objects.active = self.light

        # change location
        self.light.location = location


    def apply_global_render_settings(self, renderer='BLENDER_WORKBENCH', exposure=0.5, gamma=1.5, samples=4):

        ## some color and lighting settings
        bpy.context.scene.view_settings.view_transform = 'Standard'
        # bpy.context.space_data.shading.type = 'RENDERED'

        ## transparent background
        bpy.context.scene.render.film_transparent = True
        bpy.context.scene.render.image_settings.color_mode = 'RGBA'

        ## use GPU for rendering
        bpy.context.scene.render.engine = renderer
        if renderer == 'CYCLES':
            bpy.context.preferences.addons["cycles"].preferences.compute_device_type = "CUDA"  # or "OPENCL"
            bpy.context.scene.cycles.device = "GPU"
            ## basically the higher the sharper the render
            bpy.context.scene.cycles.samples = samples
        else:
            bpy.context.scene.cycles.device = "CPU"
            bpy.data.scenes['Scene'].display.render_aa = str(samples)



    def add_color(self,obj,remove=True):
        if remove:
            # remove color from wireframe
            col=obj.data.color_attributes.get('Col')
            obj.data.color_attributes.remove(col)
        # add new black color
        colattr = obj.data.color_attributes.new(
            name='my_color',
            type='FLOAT_COLOR',
            domain='POINT',
        )
        cols = []
        for v_index in range(len(obj.data.vertices)):
            cols += [0, 0, 0, 1]
        colattr.data.foreach_set("color", cols)

    def render_surface(self, file, out=None):

        ## get file and put it in scene collection
        bpy.ops.import_mesh.ply(filepath=file)
        obj1 = bpy.context.active_object
        if self.rotate:
            # obj.matrix_world = Matrix.Rotation(-math.pi/2,4,"Z") @ Matrix.Rotation(math.pi/2,4,"X")
            obj1.matrix_world = Matrix.Rotation(self.rotate[2], 4, "Z") \
                               @ Matrix.Rotation(self.rotate[1], 4, "Y") \
                               @ Matrix.Rotation(self.rotate[0], 4, "X")
        # remove it from scene collection
        self.scene_coll.objects.unlink(obj1)
        self.coll.objects.link(obj1)



        bpy.ops.import_mesh.ply(filepath=file)
        obj = bpy.context.active_object

        modifier = obj.modifiers.new(name='my_modifier',type='WIREFRAME')
        modifier.thickness = 0.001
        bpy.ops.object.modifier_apply(modifier='my_modifier')

        self.add_color(obj,remove=True)

        if self.rotate:
            # obj.matrix_world = Matrix.Rotation(-math.pi/2,4,"Z") @ Matrix.Rotation(math.pi/2,4,"X")
            obj.matrix_world = Matrix.Rotation(self.rotate[2], 4, "Z") \
                               @ Matrix.Rotation(self.rotate[1], 4, "Y") \
                               @ Matrix.Rotation(self.rotate[0], 4, "X")
        # remove it from scene collection
        self.scene_coll.objects.unlink(obj)
        self.coll.objects.link(obj)


        outfile = out if out else str(Path(file).with_suffix(".png"))
        bpy.context.scene.render.filepath = outfile
        bpy.ops.render.render(write_still=True)

        print("Mesh render saved to ", outfile)

        if self.remove_model:
            self.coll.objects.unlink(obj)
            self.coll.objects.unlink(obj1)



    def color_along_axis(self, points, axis=1):
        sign = np.sign(axis)
        axis = np.abs(axis) - 1
        cmap = 'jet'
        return MplColorHelper(cmap, points[:, axis].min(), points[:, axis].max()).get_rgb(sign * points[:, axis])
        # cols=MplColorHelper(cmap, accuracy.min(), accuracy.max()).get_rgb(accuracy)

    def render_point_cloud(self, file, rotate=None, out=None):


        """this is the one to use for point cloud rendering"""
        if os.path.splitext(file)[1] == ".ply":
            pcd = PyntCloud.from_file(file)
            points = pcd.points[["x", "y", "z"]].values
            normals = pcd.points[["nx", "ny", "nz"]].values
            if "red" in pcd.points.keys():
                colors = pcd.points[["red", "green", "blue"]].values
            else:
                colors = None
        elif os.path.splitext(file)[1] == ".npz":
            data = np.load(file)
            points = data["points"]
            normals = data["normals"]
            if "colors" in data.keys():
                colors = data["colors"]
            else:
                colors = None
            normals = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]
        else:
            print("ERROR: {} is not a supported file ending for point cloud rendering".format(os.path.splitext(file)[1]))

        if rotate is not None:
            points = points @ Matrix.Rotation(-rotate[0], 3, "X") @ Matrix.Rotation(-rotate[1], 3,
                                                                                         "Y") @ Matrix.Rotation(
                -rotate[2], 3, "Z")
            normals = normals @ Matrix.Rotation(-rotate[0], 3, "X") @ Matrix.Rotation(-rotate[1], 3,
                                                                                           "Y") @ Matrix.Rotation(
                -rotate[2], 3, "Z")


        # # normals to rotmat
        # marker_default_orient = [0, 0, 1]
        # angles = np.arccos(np.dot(normals, marker_default_orient))
        # cross = np.cross(normals, marker_default_orient)
        # cross = cross / np.linalg.norm(cross, axis=1)[:, np.newaxis]
        # quat = np.array([cross[:, 0], cross[:, 1], cross[:, 2], angles[:]]).transpose()
        # rots = Rotation.from_quat(quat).as_matrix()
        # # has to be rotated 180 degrees when using cone marker
        # I = np.identity(3)
        # I[0, 0] = -1
        # I[2, 2] = -1
        # rots = rots @ I
        #
        # ## colors
        colors = self.color_along_axis(points=points, axis=self.color_axis)[:, :3]
        #
        pc = "pc"
        # bplt.Scatter(points,
        #              color=colors,
        #              marker_type="cones",
        #              radius_bottom=1,
        #              radius_top=3,
        #              marker_scale=[self.marker_scale, self.marker_scale, self.marker_scale / 3],
        #              marker_rotation=rots,
        #              randomize_rotation=False,
        #              name=pc)
        bplt.Scatter(points,
                     color=colors,
                     marker_type="uv_spheres",
                     marker_scale=[self.marker_scale,self.marker_scale,self.marker_scale],
                     randomize_rotation=False,
                     name=pc)

        obj = bpy.context.scene.objects[pc]
        self.scene_coll.objects.unlink(obj)
        self.coll.objects.link(obj)

        outfile = out if out else str(Path(file).with_suffix(".png"))
        bpy.context.scene.render.filepath = outfile
        bpy.ops.render.render(write_still=True)
        print("Renderer to", outfile)

        if self.remove_model:
            self.coll.objects.unlink(obj)

    def apply_render_settings(self, mode="normal"):

        self.rotate = None

        match mode:
            case "normal":
                bpy.data.scenes['Scene'].display.shading.light = 'MATCAP'
                bpy.data.scenes['Scene'].display.shading.studio_light = 'check_normal+y.exr'
                bpy.data.scenes['Scene'].display.shading.color_type = 'OBJECT'
            case "single_color":
                bpy.data.scenes['Scene'].display.shading.light = 'STUDIO'
                bpy.data.scenes['Scene'].display.shading.studio_light = 'rim.sl'
                bpy.data.scenes['Scene'].display.shading.single_color = (0.8, 0.183968, 0)
                bpy.data.scenes['Scene'].display.shading.color_type = 'SINGLE'
            case "axis":
                bpy.data.scenes['Scene'].display.shading.light = 'STUDIO'
                bpy.data.scenes['Scene'].display.shading.studio_light = 'rim.sl'
                bpy.data.scenes['Scene'].display.shading.color_type = 'VERTEX'
            case "color":
                bpy.data.scenes['Scene'].display.shading.light = 'STUDIO'
                bpy.data.scenes['Scene'].display.shading.studio_light = 'rim.sl'
                bpy.data.scenes['Scene'].display.shading.color_type = 'VERTEX'

        self.color_axis = 2

        self.marker_scale = 0.25

        # set in blender in "output properties" -> format -> resolution
        self.resolution = (768, 1024)
        bpy.context.scene.view_settings.exposure = 0.4
        bpy.context.scene.view_settings.gamma = 1.6
        bpy.context.scene.view_settings.look = 'Medium High Contrast'



if __name__ == "__main__":

    # TODO: add a plane below the lowest point of the point cloud so I can cast a shadow on it

    remove_model = False

    path = "/home/rsulzer/cpp/psdr/example/data"
    model = "gargoyle"


    br = BlenderRender()
    br.remove_model = remove_model


    br.apply_render_settings()

    br.add_camera(os.path.join(path,model,"camera.npz"), br.resolution)

    # br.apply_global_render_settings(renderer='BLENDER_WORKBENCH', samples=5)
    br.apply_global_render_settings(renderer='CYCLES', samples=5)

    br.render_point_cloud(os.path.join(path,model,"pointcloud.ply"), rotate=[math.pi/2,0,0])

    if br.remove_model:
        bpy.data.objects.remove(br.camera, do_unlink=True)

    a=5

