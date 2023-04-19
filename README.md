# TensorAnalysis
Calculate descriptors (globularity, asphericity, acylindricity, anisotropy), from tensor data extracted using the cpptraj tool.

To use:
> python tensor.py -i <tensor.dat output>

To get the output tensor.dat file, in cpptraj:
>radgyr out tensor.dat tensor
