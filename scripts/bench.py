import matplotlib.pyplot as plt
import numpy as np
import subprocess

max_n = 6
sizes = [2**i for i in range(max_n)]
target_ts = [10 for i in range(max_n)]

def no_save():
    data = {}

    for i, s, t in zip(list(range(max_n)), sizes, target_ts):
        print(f"{s} bodies. target_t: {t}")
        out = subprocess.run(["./build/src/bench", str(s), str(t), "30"], capture_output=True)
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
        plt.plot(sizes, sub["times"] * 1e9 / sub["f_evals"], label=name)

    plt.title("accel. f. evaluation efficiency")
    plt.xlabel("bodies")
    plt.ylabel("ns / acceleration function eval")
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

    for name, sub in data.items():
        e_tot = 0
        lines = []
        with open(f"{name}.txt", "r") as f:
            lines = f.readlines()

        e_tot = float(lines[1].split(";")[-1])
        err = []
        for i in range(1, len(lines)):
            err.append((e_tot - float(lines[i].split(";")[-1])) / e_tot)

        plt.plot(range(len(err)), err, label=name)
    
    plt.title("Relative error")
    plt.xlabel("step")
    plt.ylabel("e")
    plt.grid(True)
    plt.legend()
    plt.show()

save()
        

