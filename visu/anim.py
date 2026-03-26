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
            steps = data[0]
            m1 = data[1]
            p1 = data[2:5].T
            v1 = data[5:8].T
            m2 = data[8]
            p2 = data[9:12].T
            v2 = data[12:15].T

            m_scale = max(m1[0], m2[0])

            path1 = VMobject(stroke_color=colors[i], stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(p1))
            path2 = VMobject(stroke_color=colors[i], stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(p2))
            dot1 = Dot3D(radius=DEFAULT_DOT_RADIUS * m1[0] / m_scale, color=colors[i])
            dot2 = Dot3D(radius=DEFAULT_DOT_RADIUS * m2[0] / m_scale, color=colors[i])
        
            paths.append(path1)
            paths.append(path2)
            dots.append(dot1)
            dots.append(dot2)

        # nah wtf python
        for i in range(len(dots)):
            exec(f"dots[{i}].add_updater(lambda x: x.move_to(paths[{i}].get_end()))", globals=locals(), locals=globals())

        self.add(*dots)
        self.play(
            *[Create(path, run_time=5, rate_func=linear) for path in paths]
        )
