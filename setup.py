from setuptools import setup, find_packages
setup(
    name = "MAUND",
    scripts=[
        'bin/maund.py',
    ],
    packages=find_packages(exclude=['bin']),
    install_requires=[
        "pandas",
        "python-Levenshtein",
    ],
    version = "0.3.3.1",
    description = "Miseq analysis program",
    )


