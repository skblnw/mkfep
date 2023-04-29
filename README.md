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

```
conda create -n alchemical_analysis python=2.7
conda activate alchemical_analysis
conda install -c conda-forge pymbar
conda install matplotlib
git clone https://github.com/MobleyLab/alchemical-analysis.git
cd alchemical-analysis
python setup.py install
pip install pathlib2
```

## Applying the Patch

Apply the following patch to the `alchemical_analysis.py` file:

```
> diff alchemical_analysis.py alchemical_analysis.py.backup
334c334
<          O = MBAR.computeOverlap()['matrix']
---
>          O = MBAR.computeOverlap()[2]
```


You're now ready to use the Alchemical Analysis package in your Python projects!


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed under the [LICENSE_NAME] License - see the [LICENSE](LICENSE) file for details.

## Contact

- Your Name - your.email@example.com
- Project Link: https://github.com/your-username/your-repository
