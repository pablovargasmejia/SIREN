from setuptools import setup, find_packages

setup(
    name="siren-rnai",
    version="0.1.6",
    author="Pablo Vargas Mejia",
    author_email="",
    description="SIREN: Suite for Intelligent RNAi design and Evaluation of Nucleotide sequences",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pablovargasmejia/SIREN",
    license="GPLv3",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "biopython",
        "tqdm",
        "primer3-py",
        "matplotlib",
        "numpy",
        "seaborn",       
        "scipy"
    ],
    entry_points={
    "console_scripts": [
        "SIREN = siren_rnai.siren_masterV:main",  
    ],
},
)

