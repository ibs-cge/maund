from setuptools import setup, find_packages
setup(
    name = "MAUND",
    scripts=[
        'bin/maund.py',
    ],
    install_requires=[
        "pandas",
        "python-Levenshtein",
    ],
    version = "0.3.2.1",
    description = "Miseq analysis program",
    )


