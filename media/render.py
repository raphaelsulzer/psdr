import sys

import bpy, math, os
# import pydevd_pycharm
# pydevd_pycharm.settrace('localhost', port=1090, stdoutToServer=True, stderrToServer=True)

import numpy as np
import blender_plots as bplt
# from scipy.spatial.transform import Rotation

from glob import glob
# from pyntcloud import PyntCloud

from mathutils import Vector, Matrix

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib import cm


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

    def __init__(self,settings,remove_model=False):

        self.settings = settings

        self.remove_model = remove_model

        self.rotation = None

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

    def add_moving_camera(self, cam_folder, resolution=(1024, 1024), frames=None):

        cam_files = glob(os.path.join(cam_folder,'*'))

        cam_files = sorted(cam_files)
        cam_file = os.path.join(cam_folder, cam_files[0])

        cam_data = np.load(cam_file)
        # get with C.scene.camera.location
        location = cam_data["location"]
        # get with C.scene.camera.matrix_world.to_euler()
        orientation = cam_data["orientation"]

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

        keys = np.arange(0, frames + int(frames / (len(cam_files) -1 )), int(frames / (len(cam_files) - 1)))

        for i,cf in enumerate(cam_files):

            print("Add camera {} to key {}".format(cf,keys[i]))

            cam_file = os.path.join(cam_folder, cf)
            cam_data = np.load(cam_file)
            location = cam_data["location"]
            orientation = cam_data["orientation"]


            self.camera.location = location
            self.camera.rotation_euler = orientation
            self.camera.keyframe_insert(data_path='location',frame=keys[i])
            self.camera.keyframe_insert(data_path='rotation_euler',frame=keys[i])

            # ttc=self.camera.constraints.new(type='TRACK_TO')
            # ttc.target = self.object

        bpy.context.scene.camera = self.camera


    def add_rotating_camera(self, cam_folder, resolution=(1024, 1024), frames=None):

        cam_files = glob(os.path.join(cam_folder,'*'))

        cam_files = sorted(cam_files)
        cam_file = os.path.join(cam_folder, cam_files[0])

        cam_data = np.load(cam_file)
        # get with C.scene.camera.location
        location = cam_data["location"]
        # get with C.scene.camera.matrix_world.to_euler()
        orientation = cam_data["orientation"]

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

        keys = np.arange(0, frames + int(frames / (len(cam_files) -1 )), int(frames / (len(cam_files) - 1)))

        for i,cf in enumerate(cam_files):

            print("Add camera {} to key {}".format(cf,keys[i]))

            cam_file = os.path.join(cam_folder, cf)
            cam_data = np.load(cam_file)
            location = cam_data["location"]
            orientation = cam_data["orientation"]


            self.camera.location = location
            self.camera.rotation_euler = orientation
            self.camera.keyframe_insert(data_path='location',frame=keys[i])
            self.camera.keyframe_insert(data_path='rotation_euler',frame=keys[i])

            # ttc=self.camera.constraints.new(type='TRACK_TO')
            # ttc.target = self.object

        bpy.context.scene.camera = self.camera


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


        # ttc=self.camera.constraints.new(type='TRACK_TO')
        # ttc.target = self.object

        bpy.context.scene.camera = self.camera



    def add_light(self, location, energy=1000000):

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
        # cannot get white background easily when rendering with cycles: https://www.reddit.com/r/blenderhelp/comments/azr9h0/comment/ei9lsx3/?utm_source=share&utm_medium=web3x&utm_name=web3xcss&utm_term=1&utm_content=share_button
        bpy.context.scene.render.film_transparent = True
        bpy.context.scene.render.image_settings.color_mode = 'RGBA'

        ## use GPU for rendering
        bpy.context.scene.render.engine = renderer
        if renderer == 'CYCLES':
            bpy.context.preferences.addons["cycles"].preferences.compute_device_type = "CUDA"  # or "OPENCL"
            bpy.context.scene.cycles.device = "GPU"
            # bpy.context.scene.cycles.device = "CPU"
            ## basically the higher the sharper the render
            bpy.context.scene.cycles.samples = samples
        else:
            bpy.data.scenes['Scene'].display.render_aa = str(samples)

    def add_color(self,obj,color,remove=True):
        if remove:
            # remove color from wireframe
            col=obj.data.color_attributes.get('Col')
            obj.data.color_attributes.remove(col)
        # add new black color
        colattr = obj.data.color_attributes.new(
            name='Col',
            type='FLOAT_COLOR',
            domain='POINT',
        )
        cols = []
        for v_index in range(len(obj.data.vertices)):
            cols += color # has to be with alpha value
        colattr.data.foreach_set("color", cols)


    def get_object_as_np_array(self,obj):

        # from here: https://blog.michelanders.nl/2016/02/copying-vertices-to-numpy-arrays-in_4.html

        verts = np.empty(len(obj.data.vertices) * 3, dtype=np.float64)
        obj.data.vertices.foreach_get('co', verts)

        verts=np.reshape(verts, (len(obj.data.vertices), 3))

        rot = np.array(obj.matrix_world)[:3, :3]

        return np.dot(verts,rot.transpose())


    def add_surface(self, file, rotation=None, color=None, use_vertex=False):


        ## get file and put it in scene collection
        bpy.ops.import_mesh.ply(filepath=file,use_verts=use_vertex)
        obj1 = bpy.context.active_object
        self.object = obj1


        if rotation is not None:
            obj1.matrix_world = self.rotate(None,rotation)

        if color is not None:
            self.add_color(obj1,color=color,remove=False)


        # remove it from scene collection
        self.scene_coll.objects.unlink(obj1)
        self.coll.objects.link(obj1)

        return obj1


    def add_wireframe(self, file, rotation=None):

        bpy.ops.import_mesh.ply(filepath=file)
        obj = bpy.context.active_object

        modifier = obj.modifiers.new(name='my_modifier',type='WIREFRAME')
        modifier.thickness = 0.02
        modifier.use_relative_offset = False
        modifier.use_even_offset = False
        bpy.ops.object.modifier_apply(modifier='my_modifier')

        self.add_color(obj,remove=False)

        if rotation:
            # obj.matrix_world = Matrix.Rotation(-math.pi/2,4,"Z") @ Matrix.Rotation(math.pi/2,4,"X")
            obj.matrix_world = self.rotate(obj.matrix_world,rotation)

        print("here")

        # remove it from scene collection
        self.scene_coll.objects.unlink(obj)
        self.coll.objects.link(obj)

        return obj


    def color_along_axis(self, points, axis=1):
        sign = np.sign(axis)
        axis = np.abs(axis) - 1
        cmap = 'jet'
        return MplColorHelper(cmap, points[:, axis].min(), points[:, axis].max()).get_rgb(sign * points[:, axis])
        # cols=MplColorHelper(cmap, accuracy.min(), accuracy.max()).get_rgb(accuracy)

    def rotate(self, object, rotation_vector):

        if object is not None:
            return object @ Matrix.Rotation(-rotation_vector[0], 3, "X") \
                   @ Matrix.Rotation(-rotation_vector[1], 3, "Y") \
                   @ Matrix.Rotation(   -rotation_vector[2], 3, "Z")
        else:
            return Matrix.Rotation(-rotation_vector[0], 4, "X") \
                   @ Matrix.Rotation(-rotation_vector[1], 4, "Y") \
                   @ Matrix.Rotation(-rotation_vector[2], 4, "Z")

    def add_point_cloud_bplt_scatter(self, file, rotation=None, color=None):


        """this is the one to use for point cloud rendering"""
        if os.path.splitext(file)[1] == ".ply":
            pcd = PyntCloud.from_file(file)
            points = pcd.points[["x", "y", "z"]].values
            normals = pcd.points[["nx", "ny", "nz"]].values
            normals = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]
            if "red" in pcd.points.keys():
                color = pcd.points[["red", "green", "blue"]].values

        elif os.path.splitext(file)[1] == ".npz":
            data = np.load(file)
            points = data["points"]
            normals = data["normals"]
            if "colors" in data.keys():
                color = data["colors"]
            normals = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]
        else:
            print("ERROR: {} is not a supported file ending for point cloud rendering".format(os.path.splitext(file)[1]))

        if rotation is not None:
            points = self.rotate(points,rotation)
            normals = self.rotate(normals,rotation)


        self.model_bb = np.array((points.min(axis=0),points.max(axis=0)))


        # normals to rotmat
        marker_default_orient = [0, 0, 0.99]
        angles = np.arccos(np.dot(normals, marker_default_orient))
        cross = np.cross(normals, marker_default_orient)
        # cross = cross / np.linalg.norm(cross, axis=1)[:, np.newaxis]
        quat = np.array([cross[:, 0], cross[:, 1], cross[:, 2], angles[:]]).transpose()

        rots = Rotation.from_quat(quat).as_matrix()
        # has to be rotated 180 degrees when using cone marker
        I = np.identity(3)
        I[0, 0] = -1
        I[2, 2] = -1
        rots = rots @ I
        #
        # ## colors
        if color is None:
            color = self.color_along_axis(points=points, axis=self.color_axis)[:, :3]

        #
        pc = "pc"

        print("Loaded {} points".format(points.shape[0]))

        points = points[:1000000]
        rots = rots[:1000000]

        # bplt.Scatter(points,
        #              color=color,
        #              marker_type="cones",
        #              radius_bottom=1,
        #              radius_top=3,
        #              marker_scale=[self.marker_scale, self.marker_scale, self.marker_scale / 3],
        #              # marker_scale=[self.marker_scale, self.marker_scale, self.marker_scale],
        #              marker_rotation=rots,
        #              randomize_rotation=False,
        #              name=pc)
        bplt.Scatter(points,
                     # color=color,
                     marker_type="circles",
                     marker_scale=[self.marker_scale, self.marker_scale, self.marker_scale],
                     randomize_rotation=False,
                     name=pc)


        obj = bpy.context.scene.objects[pc]
        self.object = obj
        self.scene_coll.objects.unlink(obj)
        self.coll.objects.link(obj)

        return obj



    def render(self,outfile):

        # outfile = out if out else str(Path(file).with_suffix(".png"))
        bpy.context.scene.render.filepath = outfile
        bpy.ops.render.render(write_still=True)
        print("Renderer to", outfile)


    def add_shadow_catcher(self):

        # ### cube
        # k = 5
        # scale = self.object.dimensions*k
        # bpy.ops.mesh.primitive_cube_add(size=1.0, calc_uvs=True, enter_editmode=False, align='WORLD',
        #                                 location=(0.0, 0.0, ((k-1)*self.object.dimensions[2])/2.1), rotation=(0.0, 0.0, 0.0), scale=scale)
        # plane = bpy.data.objects['Cube']

        # radius = np.max(np.abs(self.model_bb[0]-self.model_bb[1]))

        minz = np.min(self.get_object_as_np_array(self.object)[:,2])

        size = max(self.object.dimensions[0],self.object.dimensions[1])
        # location = Vector((self.object.location[0],self.object.location[1],self.object.location[2]-(self.object.dimensions[2]/2)))
        location = Vector((self.object.location[0],self.object.location[1],minz))

        bpy.ops.mesh.primitive_plane_add(size=size*4, location=location)

        plane = bpy.data.objects['Plane']
        bpy.data.objects['Plane'].is_shadow_catcher = True # this makes the plane transparent, ie ONLY a shadow catcher
        # plane.select = True
        obj = bpy.context.active_object
        self.scene_coll.objects.unlink(obj)
        self.coll.objects.link(obj)

        self.shadow_catcher = plane


    def mesh_to_points(self,object,point_size=0.0001):

        print("Apply mesh to points with point size {}".format(point_size))

        bpy.context.view_layer.objects.active = object

        # mod = object.modifiers.active
        mod = object.modifiers.new(name='my_modifier', type='NODES')

        bpy.ops.node.new_geometry_node_group_assign('INVOKE_DEFAULT')

        ntp= mod.node_group.nodes.new('GeometryNodeMeshToPoints')
        self.meshToPoints=ntp
        ntp.inputs[3].default_value = point_size

        gi = mod.node_group.nodes['Group Input']
        l = gi.outputs[0].links[0] # remove default link
        mod.node_group.links.remove(l)
        go = mod.node_group.nodes['Group Output']

        mod.node_group.links.new(gi.outputs[0], ntp.inputs[0])
        mod.node_group.links.new(ntp.outputs[0], go.inputs[0])


    def add_attribute_color_cycles(self,obj):
        # from here:
        # https://stackoverflow.com/a/69807985/20795095


        # Load material
        mymat = bpy.data.materials.get("mymat")
        if mymat is None:
            mymat = bpy.data.materials.new(name="mymat")

        mymat.use_nodes = True

        # Link material to mesh
        if obj.data.materials:
            obj.data.materials[0] = mymat
        else:
            obj.data.materials.append(mymat)

        # Get node tree from the material
        nodes = mymat.node_tree.nodes
        principled_bsdf_node = nodes.get("Principled BSDF")

        # # Get Vertex Color Node, create it if it does not exist in the current node tree
        # if not "VERTEX_COLOR" in [node.type for node in nodes]:
        #     vertex_color_node = nodes.new(type="ShaderNodeVertexColor")
        # else:
        #     vertex_color_node = nodes.get("Col")
        vertex_color_node = nodes.new(type="ShaderNodeVertexColor")

        # # Set the vertex_color layer we created at the beginning as input
        vertex_color_node.layer_name = "Col"

        # Link Vertex Color Node "Color" output to Principled BSDF Node "Base Color" input
        links = mymat.node_tree.links
        links.new(vertex_color_node.outputs[0], principled_bsdf_node.inputs[0])




    def set_point_cloud_color(self,object):

        mod = object.modifiers.active

        ntp= mod.node_group.nodes.new('GeometryNodeSetMaterial')

        gi = self.meshToPoints
        l = gi.outputs[0].links[0] # remove default link
        mod.node_group.links.remove(l)
        go = mod.node_group.nodes['Group Output']

        mod.node_group.links.new(gi.outputs[0], ntp.inputs[0])
        mod.node_group.links.new(ntp.outputs[0], go.inputs[0])

        ntp.inputs[2].default_value = bpy.data.materials.get("mymat")




    def set_freestyle_edge(self,obj,thickness=0.25,color=(0.0406086, 1, 0.038908)):
        # more or less this: https://www.youtube.com/watch?v=TE9G-Y5EQBE&ab_channel=JoshGambrell

        bpy.context.scene.render.use_freestyle = True

        # obj.select_set(True)
        bpy.context.view_layer.objects.active = obj

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')

        bpy.ops.mesh.mark_freestyle_edge(clear=False)

        bpy.data.linestyles["LineStyle"].color = color
        bpy.data.linestyles["LineStyle"].thickness = thickness

        bpy.context.view_layer.freestyle_settings.linesets["LineSet"].select_silhouette = False
        bpy.context.view_layer.freestyle_settings.linesets["LineSet"].select_crease = False
        bpy.context.view_layer.freestyle_settings.linesets["LineSet"].select_border = False
        bpy.context.view_layer.freestyle_settings.linesets["LineSet"].select_edge_mark = True



    def apply_color_settings(self, mode="normal"):


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

        self.marker_scale = 0.025

        # set in blender in "output properties" -> format -> resolution
        # self.resolution = (768, 1024)
        self.resolution = (600, 600)
        bpy.context.scene.view_settings.exposure = 0.4
        bpy.context.scene.view_settings.gamma = 1.6
        bpy.context.scene.view_settings.look = 'High Contrast'

    def add_point_cloud(self, file, rotation=None,color=None):

        # this function requires the addon "import ply as verts" from the following link to be installed
        # https://github.com/TombstoneTumbleweedArt/import-ply-as-verts

        object = self.add_surface(file, rotation=rotation, color=color,
                                use_vertex=True)
        self.mesh_to_points(object, point_size=self.settings["point_size"])
        self.add_attribute_color_cycles(object)
        self.set_point_cloud_color(object)

        return object


    def set_light_and_shadow(self,light=True,shadow_catcher=True):

        ### light and shadow
        if light:
            light = Vector((0, 0, self.object.dimensions[2] * 3.5))
            light = self.rotate(light, self.settings["light_rotation"])
            self.add_light(light, self.settings["light"])

        if shadow_catcher:
            self.add_shadow_catcher()




def get_settings_dict(model):
    model_dict = dict()
    model_dict["city"] = dict()
    model_dict["city"]["rotation"] = None
    model_dict["city"]["light"] = 20000
    model_dict["city"]["point_size"] = 0.03
    model_dict["city"]["light_rotation"] = [math.pi / 5, -math.pi / 5, 0]

    model_dict["bunny"] = dict()
    model_dict["bunny"]["camera_resolution"] = [600,600]
    model_dict["bunny"]["rotation"] = [-math.pi / 2, 0, 0]
    model_dict["bunny"]["light"] = 1
    model_dict["bunny"]["point_size"] = 0.0005
    model_dict["bunny"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]

    model_dict["anchor"] = dict()
    model_dict["anchor"]["camera_resolution"] = [400,400]
    model_dict["anchor"]["rotation"] = None
    model_dict["anchor"]["light"] = 500000
    model_dict["anchor"]["point_size"] = 0.4
    model_dict["anchor"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]

    model_dict["armadillo"] = dict()
    model_dict["armadillo"]["resolution"] = [400,400]
    model_dict["armadillo"]["rotation"] = [-math.pi / 2, 0, 0]
    model_dict["armadillo"]["light"] = 3000000
    model_dict["armadillo"]["point_size"] = 0.4
    model_dict["armadillo"]["light_rotation"] = [-math.pi / 10, -math.pi / 20, 0]

    model_dict["tarbosaurus"] = dict()
    model_dict["tarbosaurus"]["resolution"] = [400,400]
    model_dict["tarbosaurus"]["rotation"] = [math.pi, 0, 0]
    model_dict["tarbosaurus"]["light"] = 800000
    model_dict["tarbosaurus"]["point_size"] = 0.4
    model_dict["tarbosaurus"]["light_rotation"] = [-math.pi / 10, -math.pi / 20, 0]


    model_dict["default"] = dict()
    model_dict["default"]["rotation"] = None
    model_dict["default"]["light"] = 50000
    model_dict["default"]["point_size"] = 0.4
    model_dict["default"]["light_rotation"] = [math.pi / 20, math.pi / 20, 0]

    return model_dict[model]

def get_path(mode):
    if mode in ["convexes_refined", "convexes_detected",
                "pointgroups_refined", "pointgroups_detected",
                "pointcloud", "alpha",
                "pointgroups", "rectangle", "convex_hull", "alpha_shape"]:
        path = "/home/rsulzer/cpp/psdr/example/data"
    elif mode in ["colored_soup", "polygon_mesh", "dense_mesh", "polygon_mesh_detected","partition"]:
        path = "/home/rsulzer/python/compod/example/data"
    else:
        print("Mode {} does not exist".format(mode))
        raise ValueError
    return path


def render_all():
    model = "anchor"
    model = "armadillo"


    modes = ["colored_soup","polygon_mesh","dense_mesh","pointcloud","convexes_refined","convexes_refined_samples","convexes_detected","convexes_detected_samples"]
    modes = ["colored_soup","polygon_mesh","dense_mesh","pointcloud","convexes_refined","convexes_refined_samples"]
    modes = ["pointcloud","convexes_refined","convexes_refined_samples"]


    modes = ["alpha"]
    # modes = ["pointgroups","rectangle","convex_hull","alpha_shape"]

    # modes = ["convexes_detected","pointgroups_detected","convexes_refined","pointgroups_refined"]


    settings = get_settings_dict(model)

    for scale in [40,100,500,2000,10000]:
    # for scale in [""]:
        # for scale in [40]:
        for mode in modes:

            path = get_path(mode)


            frames = 150

            br = BlenderRender(settings)
            br.apply_color_settings("color")
            # br.apply_global_render_settings(renderer='BLENDER_WORKBENCH', samples=5)
            br.apply_global_render_settings(renderer='CYCLES', samples=512)


            if mode in ["colored_soup","polygon_mesh","polygon_mesh_detected","dense_mesh",
                        "convexes_refined","convexes_detected","alpha",
                        "rectangle","convex_hull","alpha_shape"]:
                object = br.add_surface(os.path.join(path,model,mode,"file{}.ply".format(scale)),color=[0.8,0.8,0.8,1],rotation=model_dict[model]["rotation"])
                br.add_attribute_color_cycles(object)
            elif mode in ["pointcloud"]:
                object = br.add_point_cloud(os.path.join(path,model,mode,"file.ply"),color=[0.8,0.8,0.8,1],rotation=model_dict[model]["rotation"])
            elif mode in ["pointgroups_refined","pointgroups_detected","pointgroups"]:
                object = br.add_point_cloud(file=os.path.join(path, model, mode, "file.ply"),rotation=model_dict[model]["rotation"])
            else:
                print("Mode {} does not exist".format(mode))
                raise ValueError



            br.set_light_and_shadow()

            ### add camera
            # br.add_moving_camera(os.path.join(path,model,"cameras"), br.resolution, frames = frames)
            br.add_camera(os.path.join(path, model, "cameras", "1.npz"), self.settings["camera_resolution"])


            if True:
                bpy.data.scenes["Scene"].frame_start = 0
                bpy.data.scenes["Scene"].frame_end = 0
                outfile = os.path.join(path,model,"renders",mode+"{}".format(scale),'i')
                bpy.context.scene.render.filepath = outfile
                bpy.ops.render.render(write_still=False,animation=True)
                print("Renderer to", outfile)

            if False:
                if mode in ["dense_mesh"]:
                    br.set_freestyle_edge(object, color=(1,0,0), thickness=0.5)
                elif mode in ["polygon_mesh","polygon_mesh_detected"]:
                    br.set_freestyle_edge(object, color=(0.0406086, 1, 0.038908), thickness=0.5)

                bpy.data.scenes["Scene"].frame_start = frames+1
                bpy.data.scenes["Scene"].frame_end = frames+1
                outfile = os.path.join(path,model,"renders",mode,'i')
                bpy.context.scene.render.filepath = outfile
                bpy.ops.render.render(write_still=False,animation=True)
                print("Renderer to", outfile)

            if False:
                bpy.data.scenes["Scene"].frame_start = frames + 1
                bpy.data.scenes["Scene"].frame_end = frames + 1
                outfile = os.path.join(path, model, "renders", mode, 'i')
                bpy.context.scene.render.filepath = outfile
                bpy.ops.render.render(write_still=False, animation=True)
                print("Renderer to", outfile)


            if br.remove_model:
                br.coll.objects.unlink(object)
                br.coll.objects.unlink(br.shadow_catcher)
                bpy.data.objects.remove(br.light, do_unlink=True)
                bpy.data.objects.remove(br.camera, do_unlink=True)

            del br


def render_bunny():

    model = "bunny"

    settings = get_settings_dict(model)

    br = BlenderRender(settings)
    br.apply_color_settings("color")
    # br.apply_global_render_settings(renderer='BLENDER_WORKBENCH', samples=5)
    br.apply_global_render_settings(renderer='CYCLES', samples=16)

    ## point cloud
    # br.add_point_cloud(os.path.join(get_path("pointcloud"), model, "pointcloud", "file.ply"),rotation=settings["rotation"])

    ## add polygons
    # object=br.add_surface(os.path.join(get_path("convexes_detected"), model, "convexes_detected", "file.ply"),rotation=settings["rotation"])
    # br.add_attribute_color_cycles(object)

    ## add partition step by step
    # for i in range()
    # object=br.add_surface(os.path.join(get_path("convexes_detected"), model, "convexes_detected", "file.ply"),rotation=settings["rotation"])
    # br.add_attribute_color_cycles(object)

    ## need to probably use wireframe modifier again for displaying the partition
    object=br.add_surface(os.path.join(get_path("partition"), model, "partition", "file.ply"),rotation=settings["rotation"])
    br.add_attribute_color_cycles(object)

    br.set_freestyle_edge(object, color=(1, 0, 0), thickness=2)

    br.set_light_and_shadow()


    ## add camera
    # br.add_moving_camera(os.path.join(path,model,"cameras"), br.resolution, frames = frames)
    br.add_camera(os.path.join(get_path("pointcloud"), model, "cameras", "1.npz"), settings["camera_resolution"])


if __name__ == "__main__":

    # render_all()
    render_bunny()












