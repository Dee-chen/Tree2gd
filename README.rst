.. image:: https://img.shields.io/pypi/v/Tree2gd.svg
   :alt: Tree2gd on the Python Package Index (PyPI)
   :target: https://pypi.python.org/pypi/Tree2gd


Tree2gd
=====================
Tree2GD provides an integrated pipeline to identify WGD events, with friendful commands in one-step or multiple steps, with smart quality control in custom dataset, with multithreading design costing low time, with well performance in detect WGD signals, and with advanced visualization of GDs and Ks peaks.


Python Requirements
===================
We currently recommend using Python 3.8 from http://www.python.org
Tree2GD is currently supported and tested on the following Python
implementations:

- Python 3.6, 3.7, 3.8 -- see http://www.python.org

- Pip3.6>=v19.2.3 -- see https://pypi.org/project/pip/
Installing by ``sudo apt-get install python3-pip`` or with get-pip.py using::

   curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
   python3 get-pip.py

R Requirements
===================
In order to complete the final drawing result display，We currently recommend using R>=v4.0.0
from https://www.r-project.org

Installation
===================
-Installation From Pypi
You can quickly install Tree2gd through the following command, and automatically install the python packages it depends on::
    pip3 install Tree2gd [--user] #You may need the --user parameter if you do not have administrator rights
If you need to use a specific python3 path to install and use Tree2gd, you can replace the above ``pip3`` with ``/THE/PATH/OF/YOUR/PYTHON3 -m pip``
-Installation From Source
You can download and decompress our source code, or fetch it using git.
Now change directory to the Tree2gd source code folder and run::
    python3 setup.py build
    python3 setup.py test
    python3 setup.py install [--user] #You may need the --user parameter if you do not have administrator rights
