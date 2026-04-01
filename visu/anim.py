from manim import *
import numpy as np

colors = [BLUE, GOLD, TEAL, RED, GREEN, MAROON, YELLOW, PURPLE, WHITE]

files = []

with open("../index.txt", "r") as f:
    for l in f.read().strip().split("\n"):
        files.append(l)

class Haha(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(
            (-4, 4, 1),
            (-4, 4, 1),
            (-4, 4, 1),
        )
        self.add(axes)
        self.set_camera_orientation(phi=30*DEGREES, theta=-45*DEGREES)

        texts = []
        legends = []
        paths = []
        dots = []

        for i, f in enumerate(files):
            color = colors[i]

            square = Rectangle(color, 0.25, 0.5).set_fill(color, 1)
            if not len(legends):
                square.to_corner(UL)
            else:
                square.move_to(legends[i-1].get_bottom() + DOWN * 0.25)
            legends.append(square)

            text = Text(f, font_size=DEFAULT_FONT_SIZE * 0.25, should_center=False)
            text.move_to(square.get_right() + RIGHT * 0.8)
            texts.append(text)

            data = np.loadtxt(f"../{f}.txt", delimiter=";", skiprows=1).T
            bodies = int(data[0][0])
            steps = data[1]
            offset = 2
            m_ratio = 5
            for j in range(bodies):
                m = data[offset+0]
                x = data[offset+1:offset+4].T
                v = data[offset+4:offset+7].T
                offset += 7

                path = VMobject(stroke_color=color, stroke_opacity=0.7, stroke_width=2).set_points_as_corners(axes.c2p(x))
                dot = Dot3D(radius=DEFAULT_DOT_RADIUS * m[0] / m_ratio, color=color)
            
                paths.append(path)
                dots.append(dot)

        # nah wtf python
        for i in range(len(dots)):
            exec(f"dots[{i}].add_updater(lambda x: x.move_to(paths[{i}].get_end()))", locals(), globals())

        self.add_fixed_in_frame_mobjects(*legends)
        self.add_fixed_in_frame_mobjects(*texts)
        self.add(*dots)
        self.begin_ambient_camera_rotation(rate=0.5)
    #    self.begin_3dillusion_camera_rotation(rate=2)
        self.play(
            *[Create(path, run_time=5, rate_func=linear) for path in paths],
        )
