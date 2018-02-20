Motility Feature Quantification
===============================

Heteromotility feature quantification can be performed on the full length of a
cell's path, or along subpaths. This tutorial document describes how to quantify
features across the entire track length, ``heteromotility``'s default behavior.

Simulated random walks and Levy fliers have been provided in the ``demo/`` directory
of the `Github page <https://github.com/cellgeometry/heteromotility>`_.

Input Data Format
-----------------

Track input formats may take one of the following two forms:

**Tracks X and Y CSV**

Heteromotility accepts object paths as ``N x T`` CSVs, where ``N`` is the number
of samples and ``T`` is the number of time steps. One CSV encodes an object's
position in the ``X`` dimension, and another in the ``Y`` dimension.

Each row describes a single sample, and each column contains the objects
position in the ``X`` or ``Y`` dimension at the corresponding time step, ordered ``0``
to ``T``, left to right.

For instance, ``tracksX.csv`` may contain the following:

.. code-block:: bash

    0, 5, 4, 5, 7, 5, 10, ...
    38, 51, 42, 38, 41, 43, ...

and likewise for ``tracksY.csv``.

**Pickled Cell Paths Object**

Heteromotility also accepts a pickled Python dictionary of object paths as an
input. The dictionary should be keyed by a unique object identifier (i.e.
numbers, names), with each key corresponding to a sequential list of XY-point

.. code-block:: python

  object_paths = {
                  obj1 : [(x1,y1), (x2,y2), (x3,y3)...],
                  obj2 : [(x1,y1), (x2,y2), (x3,y3)...],
                  ...
                  }

Feature Extraction
------------------

Simulation data provided in ``demo/`` are formated as CSVs, as outlined above.

To extract features, simply call ``heteromotility`` from the command line, specifying an
``output_path``, as well as locations for ``--tracksX`` and ``tracksY``.

.. code-block:: bash

  heteromotility demo/ --tracksX demo/rw_x.csv --tracksY demo/rw_y.csv

By default, ``heteromotility`` will export a CSV named ``motility_statistics.csv`` to the specified
``output_path`` (the first CLI argument). The output name can be altered with the ``--output_suffix`` flag.

Output Data Format
------------------

The output data is an ``N x M+2`` matrix, where ``N`` is the number of input
paths and ``M`` is the number of ``heteromotility`` features, currently ``M =
79`` using default settings. There are two ``id`` variables in the first two
columns of the matrix: ``Well/XY`` specifying the output directory, and
``cell_id`` specifying the internal unique ``id`` provided to a particular
track.

In static analysis mode, ``cell_id`` is simply an integer, ``[0, N]``.
``Well/XY`` is useful for determining cell provenance when multiple output CSVs
have been concatenated together.
