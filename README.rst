.. image:: https://img.shields.io/pypi/v/Tree2gd.svg
   :alt: Tree2gd on the Python Package Index (PyPI)
   :target: https://pypi.python.org/pypi/Tree2gd


Tree2gd
=====================
Tree2GD provides an integrated pipeline to identify WGD events, with friendful commands in one-step or multiple steps,
with smart quality control in custom dataset, with multithreading design costing low time, with well performance in detect WGD signals,
 and with advanced visualization of GDs and Ks peaks.


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
- Installation From Pypi
Y
ou can quickly install Tree2gd through the following command, and automatically install the python packages it depends on::
    pip3 install Tree2gd [--user] #You may need the --user parameter if you do not have administrator rights
If you need to use a specific python3 path to install and use Tree2gd, you can replace the above ``pip3`` with ``/THE/PATH/OF/YOUR/PYTHON3 -m pip``

- Installation From Source
You can download and decompress our source code, or fetch it using git.
Now change directory to the Tree2gd source code folder and run::
    python3 setup.py build
    python3 setup.py test
    python3 setup.py install [--user] #You may need the --user parameter if you do not have administrator rights


Testing
===================
After the installation is successful, the main program command ``Tree2gd`` and test command ``Tree2gd_test`` will be added to your system.
You can first check whether Tree2gd is installed successfully by running the following command::
    $ Tree2gd -h
If the system feeds back its corresponding parameter description, congratulations on the correct installation of Tree2gd to your system!
Next, we **strongly recommend** that you run ``Tree2gd_test`` to use the data we have prepared for a quick and complete Tree2gd test, because you can get the following benefits:
    - Check whether the pre-compiled version of the software we use by default is suitable for your system, and replace the unavailable ones with the configuration file (see the following instructions for configuration).
    - When using for the first time, in the final drawing part, we will spend a few minutes to install several dependent packages in R. After that, the formal use will be faster and more convenient.
    - After the user modifies the configuration file, you can add own new settings through the ``-config`` parameter of the command to test, and quickly detect that the new configuration can run successfully.
