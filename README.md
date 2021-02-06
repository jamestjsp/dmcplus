# DMCPplus
This is a Python project to parse Aspen DMCPlus model file as a Python dictionary and optionally dump as a json file. Some helper functions also provided to work with Finite step response (FSR) models extracted from model file. I have a plan to add more helper functions to convert MDL file to transfer functions, state space model and so on.

## MDL file

MDL file is a delimited text file with a predefined data structure. The content of the file can be summarized as follows.

1. mdl file comments
2. Number of files and the file names
3. Number of independents
4. Number of dependents
5. Number of coefficients
6. Time to steady-state
7. DMCplus model style flag
8. Tagnames
9. is Ramp flag
10. DC gain
11. Finite step response (FSR) cause and effect curves.

## Project dependencies
To parse mdl you may use mdl.py and has no dependency on the following packages but to use utility functions you may need the following libraries.
1. Numpy
2. Scipy
3. Matplotlib

## Copyright
Aspen DMCplus is registered trademarks of Aspen Technology, Inc., Cambridge, MA.
