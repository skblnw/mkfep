# mkfep

Automatic FEP setup for Gromacs

## Usage

### Preparing the environment

Make sure you have the following files and directories in your working directory:

- `md.gro` or `md.pdb`: Structure file of your system
- `pdb2fep.py`: Behind-the-scenes script to generate GROMACS topology files
- `run`: Script to execute the simulation, largely depends on user's local envrionment
- `mdp`: Directory containing MDP files for energy minimization, NVT equilibration, and MD production

### Running the script

Execute the script as follows:

```bash
# always make sure you activate your pmx environment
conda activate pmx
# look at your GMXLIB to make sure you have access to mutated ff
echo $GMXLIB

# execute the main script
python mkfep_prepare.py
```

### Analyzing the results

We have prepared an analysis script (alchemical_analysis_script_py27.py) which utilizes alchemical-analysis developed by Mobley Lab (https://github.com/MobleyLab/alchemical-analysis). 

To do: we are in the process of migrating to alchemlyb (https://github.com/alchemistry/alchemlyb). 

## Package Installation

### Steps

Download miniconda or anaconda (e.g. Miniconda3-py38_23.1.0-1-Linux-x86_64.sh)
```bash
bash Miniconda3-py38_23.1.0-1-Linux-x86_64.sh
```

#### Install pmx
```bash
conda create -n pmx python=3.8 numpy scipy matplotlib
conda activate pmx

git clone https://github.com/deGrootLab/pmx pmx
cd pmx
git checkout develop
python setup.py install
```

#### Install alchemical_analysis

***Python 2 is required***

We strongly suggest that you create a new conda environment for alchemical-analysis as it requires Python 2

```bash
conda create -n alchemical_analysis python=2.7
conda activate alchemical_analysis

conda install -c conda-forge pymbar
conda install matplotlib
git clone https://github.com/MobleyLab/alchemical-analysis.git
cd alchemical-analysis
python setup.py install
pip install pathlib2
```

### Applying the Patch

Apply the following patch to the `alchemical_analysis.py` file:

```bash
> diff alchemical_analysis.py alchemical_analysis.py.backup
334c334
<          O = MBAR.computeOverlap()['matrix']
---
>          O = MBAR.computeOverlap()[2]
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed under the [LICENSE_NAME] License - see the [LICENSE](LICENSE) file for details.

## Contact

- Kev - kevin@skblnw.com
- Project Link: https://github.com/skblnw/mkfep
