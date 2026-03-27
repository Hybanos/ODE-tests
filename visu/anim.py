from manim import *
import numpy as np

files = []

with open("../index.txt", "r") as f:
    for l in f.read().strip().split("\n"):
        files.append(f"../{l}.txt")

class Haha(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(
            (-4, 4, 1),
            (-4, 4, 1),
            (-4, 4, 1),
        )
        self.add(axes)
        self.set_camera_orientation(phi=30*DEGREES, theta=-45*DEGREES)

        paths = []
        dots = []
        colors = color_gradient([BLUE, ORANGE], len(files)) 

        for i, f in enumerate(files):
            data = np.loadtxt(f, delimiter=";", skiprows=1).T
            bodies = int(data[0][0])
            steps = data[1]
            offset = 2
            m_ratio = 5
            for j in range(bodies):
                m = data[offset+0]
                x = data[offset+1:offset+4].T
                v = data[offset+4:offset+7].T
                offset += 7

                path = VMobject(stroke_color=colors[i], stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(x))
                dot = Dot3D(radius=DEFAULT_DOT_RADIUS * m[0] / m_ratio, color=colors[i])
            
                paths.append(path)
                dots.append(dot)

        # nah wtf python
        for i in range(len(dots)):
            exec(f"dots[{i}].add_updater(lambda x: x.move_to(paths[{i}].get_end()))", globals=locals(), locals=globals())

        self.add(*dots)
        self.begin_ambient_camera_rotation(rate=0.5)
    #    self.begin_3dillusion_camera_rotation(rate=2)
        self.play(
            *[Create(path, run_time=5, rate_func=linear) for path in paths],
        )
