import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import subprocess

max_n = 6
sizes = [2**i for i in range(max_n)]
target_t = 10

colors = mcolors.TABLEAU_COLORS.values()

def no_save():
    data = {}

    for i, s in zip(list(range(max_n)), sizes):
        print(f"{s} bodies. target_t: {target_t}")
        out = subprocess.run(["./build/src/bench", str(s), str(target_t), "30"], capture_output=True)
        strout = out.stdout.decode()
        for l in strout.strip().split("\n"):
            print(l)
            name, time, steps, f_evals = l.split(";")
            time = float(time.strip()[:-1])
            steps = int(steps)
            f_evals = int(f_evals)
            if name not in data:
                data[name] = {
                    "times": np.zeros(max_n),
                    "steps": np.zeros(max_n),
                    "f_evals": np.zeros(max_n),
                }
            data[name]["times"][i] = time
            data[name]["steps"][i] = steps
            data[name]["f_evals"][i] = f_evals
        print("")

    for name, sub in data.items():
        plt.plot(sizes, sub["times"] * 1e9 / target_t, label=name)

    plt.title("accel. f. evaluation efficiency")
    plt.yscale("log", base=10)
    plt.xlabel("bodies")
    plt.ylabel("ns / final t")
    plt.grid(True)
    plt.xscale("log", base=2)
    plt.legend()
    plt.show()

def save():
    data = {}

    sizes = [10]
    target_ts = [10]

    for i, s, t in zip(list(range(max_n)), sizes, target_ts):
        print(f"{s} bodies. target_t: {t}")
        out = subprocess.run(["./build/src/bench", str(s), str(t), "30", "true"], capture_output=True)
        strout = out.stdout.decode()
        for l in strout.strip().split("\n"):
            print(l)
            name, time, steps, f_evals = l.split(";")
            time = float(time.strip()[:-1])
            steps = int(steps)
            f_evals = int(f_evals)
            if name not in data:
                data[name] = {
                    "times": np.zeros(max_n),
                    "steps": np.zeros(max_n),
                    "f_evals": np.zeros(max_n),
                }
            data[name]["times"][i] = time
            data[name]["steps"][i] = steps
            data[name]["f_evals"][i] = f_evals
        print("")

    for name, color in zip(data.keys(), colors):
        lines = []
        with open(f"{name}.txt", "r") as f:
            lines = f.readlines()

        e_energy = float(lines[1].split(";")[-1])
        e_momentum = float(lines[1].split(";")[-4])
        err_energy = []
        err_momentum = []
        x = []
        for i in range(1, len(lines)):
            err_energy.append((float(lines[i].split(";")[-1]) - e_energy) / abs(e_energy))
            err_momentum.append((float(lines[i].split(";")[-4]) - e_momentum) / abs(e_momentum))
            x.append(float(lines[i].split(";")[2]))

        plt.step(x, err_energy, label=name, color=color)
        plt.step(x, err_momentum, label=f"{name} - angular mom.", linestyle="dotted", color=color)
    
    plt.title("Relative error")
    plt.xlabel("step")
    plt.ylabel("e")
    # plt.yscale("log", base=10)
    plt.grid(True)
    plt.legend()
    plt.show()

no_save()
save()