import subprocess
import os
import concurrent.futures

def run(cmds):
    subprocess.run(cmds, shell=True)
def lastFrame(directory):
    frames = [int(file.split("geo")[1].split(".")[0]) for file in os.listdir(directory) if file.startswith("geo") and file.endswith(".mat")]
    last_frame = max(frames)
    return last_frame

here = "/home/cuncheng/Dev/active_nematics_shell/matlab"
verbose = 'false'
restart = False
rank = [10]

commands = []
kappa = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
alpha = [1, 10, 100, 1000]
for r in rank:
    for k in kappa:
        for a in alpha:
            dir = here + f"/data/nematic/gut/rk{r}_k{k}_a{a}"
            start = lastFrame(dir) if restart else 0
            print(f"start: {start}") if restart else 0
            commands.append(f""" matlab -nodisplay -nosplash -nodesktop -r  "cd '{here}'; verbose = {verbose}; dir = '{dir}/'; start = {start}; p.rank = {r}; p.kappa = {k}; p.alpha = {a}; willmore_nematic ; exit" """)

max_workers = 20
with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run, cmd) for cmd in commands]
    concurrent.futures.wait(futures)