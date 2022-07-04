# n-bit reactive strategies for direct reciprocity

This repository contains files for the project ``n-bit reactive strategies
for direct reciprocity''. More specifically, it contains all the files
necessary to compile the written notes and all the files necessary for our
analysis.

### Notes

The written notes are written in Latex. To compile the document one needs
the `main.tex`, the `bibliography.bib` and the folder `static`. `static`
contains all the static files such as figures of the notes.

### Code

The analysis has been carried out using Python, MatLab and Mathematica.

The Mathematica files can all be found in `src/mathematica`. These
files were used to obtain some of the expression presented in the notes.

MatLab was used to carry out the evolutionary simulations presented in the last
section, and the Vaquero method. These files can be found in
`src/evolutionary_simulations`. The `VaqueroMethod.m` file was used for the
method, the rest of the files were used for the evolutionary simulations.

Python was used to numerically explore the Nash in the case of two-bit
strategies, for the plotting and presentation of the results. Note that
the evolutionary process has also been implemented in Python (`src/main.py`).

The analysis was done in Jupyter Notebooks. These can be found in the folder
`nbs`. One can run these notebooks as normal. No special packages was used.