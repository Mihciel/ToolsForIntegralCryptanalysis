from setuptools import find_packages, setup
setup(
    name='Tincrbell',
    packages=find_packages(include=['Tincrbell']),
    version='0.0.1',
    description='Tools for automated integral cryptanalysis',
    author='Michiel Verbauwhede',
    # license='TODO',
    install_requires=["python-sat", "sympy", "galois", "numpy", "pqdm"],
)