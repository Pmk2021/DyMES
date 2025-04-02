# DyMES

Simulation code for the pandemic model with transition function, m_eff, and r_eff functionalities.

This code differs from the regular DyMES model, as it takes the expectation over a probability distribution from 0 to m_eff, rather than 0 to N.

This code supports the brute force implementation of lambda calculation.

## Installation
To set up this repository, first clone this branch onto your local machine. Then run

```bash
pip install -e .
```

to install all dependencies

## Usage
All necessary functions are documented in the notebook `run_model.ipynb`

## TODO
1. Modify the probability distribution summation.
2. Add the lambda dynamics method.
