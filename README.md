# DyMES

Simulation code for the pandemic model with transition function, m_eff, and r_eff functionalities.

This code differs from the regular DyMES model in one important aspect: it takes the expectation over a probability distribution from 0 to `m_eff`, rather than from 0 to N. This is because, in a pandemic model, you cannot have more than `m_eff` infected individuals in a group.

This code supports the brute force implementation of lambda calculation.

## Installation
To set up this repository, first clone this branch onto your local machine. Then run

```bash
pip install -e .
```

to install all dependencies.

## Usage
All necessary functions are documented in the notebook `run_model.ipynb`

## TODO
1. Add the lambda dynamics method
