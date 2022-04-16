#!/bin/env python
import math
import vtk

class ConnectionRenderer(object):
    renderer = None
    connection_map = dict()
    process_status = dict()

    def render_force(self, process, z_scale):
        pid = process.process_GID
        if self.process_status.has_key(pid) and self.process_status[pid]['drawn']:
            return

        self.render(process, z_scale)

    def remove_force(self, process):
        pid = process.process_GID
        if self.process_status.has_key(pid) and self.process_status[pid]['drawn']:
            self.render(process, 1.0)

    def render(self, process, z_scale):
        pid = process.process_GID

        if not self.process_status.has_key(pid):
            self.process_status[pid] = {'drawn': False, 'actors': []}

        if self.process_status[pid]['drawn']:
            self.remove_connections(process)
        else:
            self.draw_connections(process, z_scale)

    def draw_connections(self, process, z_scale):
        pid = process.process_GID

        p1 = (process.center_position[0], process.center_position[1], process.spatial_extent[5] * z_scale)

        resolution = 128

        posts = self.connection_map[pid]

        actors = []
        for p in posts:
            p2 = (p.center_position[0], p.center_position[1], p.spatial_extent[5] * z_scale)
            if p1 == p2:
                continue

            arc = vtk.vtkArcSource()
            arc.SetResolution(resolution)
            arc.SetPoint1(p1)
            arc.SetPoint2(p2)

            dx = p1[0] - p2[0]
            dy = p1[1] - p2[1]
            dz = p1[2] - p2[2]
            d = math.sqrt(3.0 * (((p1[0] - p2[0]) ** 2) + ((p1[1] - p2[1]) ** 2)))
            center = ((p1[0] + p2[0]) / 2 + dx * dz / d, (p1[1] + p2[1]) / 2 + dy * dz / d, (p1[2] + p2[2]) / 2 - d / 6)
            arc.SetCenter(center)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(arc.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetLineWidth(1)
            actor.GetProperty().SetColor([1.0, 0.0, 0.0])

            actors.append(actor)


            length = 500
            w_width = length * 0.3

            d = math.sqrt(((p1[0] - p2[0]) ** 2) + ((p1[1] - p2[1]) ** 2))
            e0 = ((p2[0] - p1[0]) / (2 * d), (p2[1] - p1[1]) / (2 * d), math.sqrt(3.) / -3.)
            pw0 = (p2[0] - e0[0] * length, p2[1] - e0[1] * length, p2[2] - e0[2] * length)

            l = math.sqrt(e0[0] ** 2 + e0[1] ** 2 + (d / math.sqrt(3.)) ** 2)
            e1 = (e0[0] / l, e0[1] / l, d / (math.sqrt(3.) * l))
            e2 = (-e1[0], -e1[1], -e1[2])

            l = math.sqrt(e1[0] ** 2 + e1[1] ** 2)
            e3 = (-e1[1] / l, e1[0] / l, 0)
            e4 = (-e3[0], -e3[1], -e3[2])

            pw1 = (pw0[0] + e1[0] * w_width, pw0[1] + e1[1] * w_width, pw0[2] + e1[2] * w_width)
            pw2 = (pw0[0] + e2[0] * w_width, pw0[1] + e2[1] * w_width, pw0[2] + e2[2] * w_width)
            pw3 = (pw0[0] + e3[0] * w_width, pw0[1] + e3[1] * w_width, pw0[2] + e3[2] * w_width)
            pw4 = (pw0[0] + e4[0] * w_width, pw0[1] + e4[1] * w_width, pw0[2] + e4[2] * w_width)

            points = vtk.vtkPoints()
            points.InsertNextPoint(p2)
            points.InsertNextPoint(pw1)
            points.InsertNextPoint(pw2)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0,0)
            triangle.GetPointIds().SetId(1,1)
            triangle.GetPointIds().SetId(2,2)

            triangles = vtk.vtkCellArray()
            triangles.InsertNextCell(triangle)

            trianglePolyData = vtk.vtkPolyData()
            trianglePolyData.SetPoints(points)
            trianglePolyData.SetPolys(triangles)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(trianglePolyData)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor([1.0, 0.0, 0.0])

            actors.append(actor)


            points = vtk.vtkPoints()
            points.InsertNextPoint(p2)
            points.InsertNextPoint(pw3)
            points.InsertNextPoint(pw4)

            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0,0)
            triangle.GetPointIds().SetId(1,1)
            triangle.GetPointIds().SetId(2,2)

            triangles = vtk.vtkCellArray()
            triangles.InsertNextCell(triangle)

            trianglePolyData = vtk.vtkPolyData()
            trianglePolyData.SetPoints(points)
            trianglePolyData.SetPolys(triangles)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(trianglePolyData)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor([1.0, 0.0, 0.0])

            actors.append(actor)


        for actor in actors:
            self.renderer.AddActor(actor)

        self.process_status[pid]['actors'] = actors
        self.process_status[pid]['drawn'] = True

    def remove_connections(self, process):
        pid = process.process_GID
        for actor in self.process_status[pid]['actors']:
            self.renderer.RemoveActor(actor)
        del self.process_status[pid]

        
class UserInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    current_z_scale = 1.0
    scale_step = 0.1

    process_map = dict()
    transFilters = []
    connectionRenderer = None
    renderWindowInteractor = None

    def __init__(self, parent=None):
        self.AddObserver("RightButtonPressEvent", self.rightButtonPressEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def rightButtonPressEvent(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
        actor = picker.GetActor()

        if actor:
            pos = actor.GetPosition()
            key = (int(pos[0]), int(pos[1]))
            if self.process_map.has_key(key):
                process = self.process_map[key]
                self.connectionRenderer.render(process, self.current_z_scale)
                self.renderWindowInteractor.GetRenderWindow().Render()

    def keyPressEvent(self, obj, event):
        symbol = self.renderWindowInteractor.GetKeySym()
        if symbol == 'Up':
            self.current_z_scale += self.scale_step
        elif symbol == 'Down' and self.current_z_scale > self.scale_step:
            self.current_z_scale -= self.scale_step
        elif symbol == 'a':
            for process in self.process_map.values():
                self.connectionRenderer.render_force(process, self.current_z_scale)
            self.renderWindowInteractor.GetRenderWindow().Render()
        elif symbol == 'z':
            for process in self.process_map.values():
                self.connectionRenderer.remove_force(process)
            self.renderWindowInteractor.GetRenderWindow().Render()
        elif symbol == 'h':
            print "###Tips of network setting visualization###"
            print "Up: expand z scale of the network"
            print "Down: reduce z scale of the network"
            print "Mouse left bottom drag: rotate the network"
            print "Mouse center bottom drag: shit the network"
            print "Mouse right bottom click: show arrows cliked tile"
            print "a: show all communication links as arrows"
            print "z: remove arrows of all communication links"
            print "q: quit visualization"
            print "h: Show Tips"

        for transFilter in self.transFilters:
            transform = vtk.vtkTransform()
            transform.Scale(1., 1., self.current_z_scale)
            transFilter.SetTransform(transform)

        self.renderWindowInteractor.GetRenderWindow().Render()

def visualize_network(process, sd):

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    userInteractorStyle = UserInteractorStyle()
    connectionRenderer = ConnectionRenderer()

    userInteractorStyle.renderWindowInteractor = renderWindowInteractor
    userInteractorStyle.connectionRenderer = connectionRenderer

    connectionRenderer.renderer = renderer

    process_hash = dict()
    process_map_from_position = dict()
    for p in process:
        process_hash[p.process_GID] = p
        process_map_from_position[(int(p.center_position[0]), int(p.center_position[1]))] = p
    userInteractorStyle.process_map = process_map_from_position

    connection_map = dict()
    for p in process:
        for pre in p.intra_regional_connection_process:
            if not connection_map.has_key(pre):
                connection_map[pre] = []
            connection_map[pre].append(p)
        for pre in p.inter_regional_connection_pre_process:
            if not connection_map.has_key(pre):
                connection_map[pre] = []
            connection_map[pre].append(p)
        for post in p.inter_regional_connection_post_process:
            if not connection_map.has_key(p.process_GID):
                connection_map[p.process_GID] = []
            connection_map[p.process_GID].append(process_hash[post])
    connectionRenderer.connection_map = connection_map

    opacity = 0.9
    layer_colors = [[0.012, 0.663, 0.957], [0.000, 0.737, 0.831], [0.000, 0.588, 0.533], [0.298, 0.686, 0.314], [0.545, 0.765, 0.290], [0.804, 0.863, 0.224], [1.000, 0.922, 0.231]]
    region_colors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]

    actors = []
    transFilters = []
    regions = []

    for p in process:
        if p.region_name not in regions:
            regions.append(p.region_name)
        region_color = region_colors[regions.index(p.region_name)]

        x0 = p.spatial_extent[0]
        y0 = p.spatial_extent[1]
        z0 = p.spatial_extent[2]
        x1 = p.spatial_extent[3]
        y1 = p.spatial_extent[4]
        z1 = p.spatial_extent[5]
        width = x1 - x0
        depth = y1 - y0
        hw = width / 2
        hd = depth / 2


        for i in range(8):
            xa = x0 if (i+1) % 4 < 2 else x1
            xb = x0 if i % 4 >= 2 else x1
            ya = y0 if i % 4 < 2 else y1
            yb = y0 if (i+1) % 4 < 2 else y1
            za = z0 if i / 4 < 1 else z1
            zb = z0 if i / 4 < 1 else z1

            line = vtk.vtkLineSource()
            line.SetPoint1(xa, ya, za)
            line.SetPoint2(xb, yb, zb)

            transform = vtk.vtkTransform()
            transform.Scale(1., 1., 1.)

            transFilter = vtk.vtkTransformFilter()
            transFilter.SetInputConnection(line.GetOutputPort())
            transFilter.SetTransform(transform)

            appendPolyData = vtk.vtkAppendPolyData()
            appendPolyData.AddInputConnection(transFilter.GetOutputPort())

            depthSortPolyData = vtk.vtkDepthSortPolyData()
            depthSortPolyData.SetInputConnection(appendPolyData.GetOutputPort())
            depthSortPolyData.SetCamera(renderer.GetActiveCamera())

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(depthSortPolyData.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetLineWidth(4)
            actor.GetProperty().SetColor(region_color)
            actors.append(actor)
            transFilters.append(transFilter)

        for i, extent in enumerate(p.subregion_positions):
            z0 = extent[2]
            z1 = extent[5]

            appendPolyData = vtk.vtkAppendPolyData()

            vertices_list = [
                [(-hw, -hd, z0), (hw, -hd, z0), (hw, hd, z0), (-hw, hd, z0)],
                [(-hw, -hd, z0), (hw, -hd, z0), (hw, -hd, z1), (-hw, -hd, z1)],
                [(hw, -hd, z0), (hw, hd, z0), (hw, hd, z1), (hw, -hd, z1)],
                [(hw, hd, z0), (-hw, hd, z0), (-hw, hd, z1), (hw, hd, z1)],
                [(-hw, hd, z0), (-hw, -hd, z0), (-hw, -hd, z1), (-hw, hd, z1)],
                [(-hw, -hd, z1), (hw, -hd, z1), (hw, hd, z1), (-hw, hd, z1)]]

            for vertices in vertices_list:
                points = vtk.vtkPoints()
                points.InsertNextPoint(vertices[0])
                points.InsertNextPoint(vertices[1])
                points.InsertNextPoint(vertices[2])
                points.InsertNextPoint(vertices[3])

                square = vtk.vtkPolygon()
                square.GetPointIds().SetNumberOfIds(4)
                square.GetPointIds().SetId(0,0)
                square.GetPointIds().SetId(1,1)
                square.GetPointIds().SetId(2,2)
                square.GetPointIds().SetId(3,3)

                squares = vtk.vtkCellArray()
                squares.InsertNextCell(square)

                polyData = vtk.vtkPolyData()
                polyData.SetPoints(points)
                polyData.SetPolys(squares)

                appendPolyData.AddInputData(polyData)

            depthSortPolyData = vtk.vtkDepthSortPolyData()
            depthSortPolyData.SetInputConnection(appendPolyData.GetOutputPort())
            depthSortPolyData.SetCamera(renderer.GetActiveCamera())

            transform = vtk.vtkTransform()
            transform.Scale(1., 1., 1.)

            transFilter = vtk.vtkTransformFilter()
            transFilter.SetInputConnection(depthSortPolyData.GetOutputPort())
            transFilter.SetTransform(transform)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transFilter.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(layer_colors[i])
            actor.GetProperty().SetOpacity(opacity)
            actor.SetPosition((x0 + x1) / 2, (y0 + y1) / 2, 0)

            actors.append(actor)
            transFilters.append(transFilter)

    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1280, 960)
    renderWindowInteractor.SetRenderWindow(renderWindow)

    userInteractorStyle.SetDefaultRenderer(renderer)
    userInteractorStyle.transFilters = transFilters
    renderWindowInteractor.SetInteractorStyle(userInteractorStyle)

    for actor in actors:
        renderer.AddActor(actor)

    renderer.SetBackground([1.0, 1.0, 1.0])
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(0)
    renderer.GetActiveCamera().Elevation(330)

    renderWindow.Render()

    renderWindowInteractor.Start()
