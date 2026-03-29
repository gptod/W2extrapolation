import sys
import shutil
import subprocess

# Pick an available Fortran compiler
compiler = next((c for c in ["gfortran", "ifort", "flang"] if shutil.which(c)), None)
if compiler is None:
    raise RuntimeError("No Fortran compiler found.")

# Compile modules
subprocess.run(
    f"{sys.executable} -m numpy.f2py --fcompiler={compiler} -c global.f90 ot_sinkhorn.f90 -m ot_sinkhorn",
    shell=True,
    check=True
)

subprocess.run(
    f"{sys.executable} -m numpy.f2py --fcompiler={compiler} -c global.f90 ot_sinkhorn.f90 w2extrapolation.f90 -m w2extrapolation",
    shell=True,
    check=True
)

print("✅ Compilation done")

