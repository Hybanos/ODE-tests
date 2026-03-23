from manim import *
import numpy as np

data = np.loadtxt("../out.txt", delimiter=";").T

steps = data[0]
m1 = data[1]
p1 = data[2:5].T
v1 = data[5:8].T
m2 = data[8]
p2 = data[9:12].T
v2 = data[12:15].T

class Haha(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(
            (-500, 500, 100),
            (-500, 500, 100),
            (-500, 500, 100),
        )
        self.add(axes)

        self.set_camera_orientation(phi=30*DEGREES, theta=-45*DEGREES)

        path1 = VMobject(stroke_color=BLUE, stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(p1))
        path2 = VMobject(stroke_color=ORANGE, stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(p2))
        dot1 = Dot3D(radius=DEFAULT_DOT_RADIUS * m1[0] / 2, color=BLUE)
        dot2 = Dot3D(radius=DEFAULT_DOT_RADIUS * m2[0] / 2, color=ORANGE)
    
        dot1.add_updater(lambda x: x.move_to(path1.get_end()))
        dot2.add_updater(lambda x: x.move_to(path2.get_end()))

        self.add(dot1, dot2)
        self.play(Create(path1, run_time=10, rate_func=linear), Create(path2, run_time=10, rate_func=linear))
        # self.play(Create(path1, runtime=10))
