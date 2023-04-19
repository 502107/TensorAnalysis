# TensorAnalysis
Calculate descriptors (globularity, asphericity, acylindricity, anisotropy), from gyration tensor data extracted using the cpptraj tool. Not the most efficient way, but gets the job done.

To use:
> python tensor.py -i <tensor.dat output>

To get the output tensor.dat file, in cpptraj:
>radgyr out tensor.dat tensor
