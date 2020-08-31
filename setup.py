#!/usr/bin/python3.5

try:
    from setuptools import setup
    from setuptools import Command
    from setuptools import Extension
except ImportError:
    sys.exit(
        "We need the Python library setuptools to be installed. "
        "Try runnning: python -m ensurepip"
    )

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='Tree2gd',
    version='1.0.31',
    packages=['tree2gd','software','example_data'],
    license='GPL',
    long_description=long_description,
    author='Duoyuan chen',
    author_email='371761145@qq.com',
    description='Tree2gd',
    py_modules=['tree2gd_main'],
    include_package_data=True,
    url='https://github.com/Dee-chen/Fa2gd',
    install_requires=[
        'configparser',
        'matplotlib',
        'argparse',
        'coloredlogs',
        'pyecharts',
        'bio'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.6',
    entry_points='''
        [console_scripts]
        Tree2gd=tree2gd_main:main
        Tree2gd_test=tree2gd_main:test
    ''',
)
