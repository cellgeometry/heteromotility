Installation
============

Heteromotility depends on ``>=python3.6``. The tool will operate using ``>=python2.7``, but this is not supported.
We recommend the Anaconda Python distribution.

Heteromotility can be installed from the Python Package Index with pip

.. code-block:: bash

  pip install heteromotility

or

Simply clone this repository and run ``setup.py`` in the standard manner.

.. code-block:: bash

  git clone https://github.com/cellgeometry/heteromotility
  cd heteromotility
  python setup.py install
  heteromotility -h
  usage: Calculate motility features from cell locations or cell paths
  [-h] [--seg SEG] [--exttrack EXTTRACK] input_dir output_dir

Both methods of package installation will add an alias ``heteromotility`` to your ``PATH`` environment variable.

As with all research software, we recommend installing Heteromotility inside of a virtual environment.

.. code-block:: bash

  virtualenv heteromotility_env/
  source heteromotility_env/bin/activate
  pip install heteromotility
