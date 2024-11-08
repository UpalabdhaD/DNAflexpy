from setuptools import setup, find_packages

setup(
    name="DNAflexpy",
    version="0.1",
    description="A tool for trinucleotide feature profiling in DNA sequences",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.0",
        "pyyaml>=5.4"
    ],
    entry_points={
        "console_scripts": [
            "DNAflexpy=DNAflexpy.cli:main"  # Replace core:main with the actual module and function you want to execute
        ]
    },
)
