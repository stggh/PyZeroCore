from setuptools import setup, find_packages

setup(
    name = "PyZeroCoreContribution",
    version = "0.1",
    description='A Python package to evaluate photodetachment cross sections and anisotropy parameters using the zero-core contribution algorithms of Stehman and Woo',
    packages=find_packages(),
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Intended Audience :: Developers',
      'Topic :: Scientific/Engineering :: Mathematics',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Scientific/Engineering :: Chemistry',
      'Topic :: Software Development :: Libraries :: Python Modules',
      'Programming Language :: Python :: 3.6',
      ],
)

