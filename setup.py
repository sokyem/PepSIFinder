from setuptools import setup, find_packages

setup(
    name="PepSIFinder",
    version="0.1.0",
    description="A package for identifying peptide isomer candidates based retention time(multiple peaks) and ion mobility data (difference between rt peaks)",
    author="Samuel Okyem",
    author_email="okyemsamuel@gmail.com",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.0.0",
        "numpy>=1.18.0"
    ],
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    include_package_data=True,
)
