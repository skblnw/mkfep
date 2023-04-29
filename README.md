# mkfep

Automatic FEP setup for Gromacs

## Features

- Feature 1
- Feature 2
- Feature 3

## Usage

### Preparing the environment

Make sure you have the following files and directories in your working directory:

- `md.gro` or `md.pdb`: Structure file of your protein
- `mkgmx_pdb2fep.py`: Script to generate GROMACS topology files
- `run`: Script to execute the simulation
- `mdp`: Directory containing MDP files for energy minimization, NVT equilibration, and MD production

### Running the script

Execute the script as follows:

```bash
python main.py
```

# Alchemical Analysis Package Installation

This README provides the installation steps for the Alchemical Analysis Python package.

## Installation Steps

1. Create a new conda environment with Python 2.7:

```
conda create -n alchemical_analysis python=2.7
```

2. Activate the newly created conda environment:

```
conda activate alchemical_analysis
```

3. Install pymbar from the conda-forge channel:

```
conda install -c conda-forge pymbar
```

4. Install matplotlib:

```
conda install matplotlib
```

5. Clone the alchemical-analysis repository:

```
git clone https://github.com/MobleyLab/alchemical-analysis.git
```

6. Change to the cloned directory:

```
cd alchemical-analysis
```

7. Install the package:

```
python setup.py install
```

8. Install pathlib2:

```
pip install pathlib2
```

## Applying the Patch

Apply the following patch to the `alchemical_analysis.py` file:

❯ diff alchemical_analysis.py alchemical_analysis.py.backup
334c334
<          O = MBAR.computeOverlap()['matrix']
---
>          O = MBAR.computeOverlap()[2]


You're now ready to use the Alchemical Analysis package in your Python projects!


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed under the [LICENSE_NAME] License - see the [LICENSE](LICENSE) file for details.

## Contact

- Your Name - your.email@example.com
- Project Link: https://github.com/your-username/your-repository
