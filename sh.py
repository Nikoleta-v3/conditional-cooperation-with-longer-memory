import numpy as np

for c in np.linspace(0, 1, 11):
    string = f"""#!/bin/bash
#SBATCH -p amd  # partition (queue)
#SBATCH -J c_{c} # job name
#SBATCH -N 1 # number of nodes, use 1-1 if you need exactly one node
#SBATCH -n 1 # number of cores
#SBATCH -t 1-00:00  # time (D-HH:MM)
#SBATCH -o c_{round(c, 2)}.out # STDOUT
#SBATCH -e c_{round(c, 2)}.err # STDERR

module load python/3.6.6
source ~/reactive/bin/activate

python main.py {c}
    """

    with open(f"sh/c_{c}.sh", "w") as text_file:
        text_file.write(string)