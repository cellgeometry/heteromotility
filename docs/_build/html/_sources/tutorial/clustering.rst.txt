Unsupervised Clustering
=======================

Heteromotility features can be used to identify heterogeneous **motility
states** in a cell population by unsupervised clustering. Tools to perform
unsupervised clustering are available in the ``analysis/`` directory of the
`Github page <https://github.com/cellgeometry/heteromotility>`_. Here, we
demonstrate two methods of clustering cells based on ``heteromotility`` feature
outputs in ``R``.

To begin, we generate ``heteromotility`` feature data from simulated random walks and Levy fliers.

.. code-block:: bash

  heteromotility demo/ --tracksX demo/rw_x.csv --tracksY demo/rw_y.csv --output_suffix rw
  heteromotility demo/ --tracksX demo/pf_x.csv --tracksY demo/pf_y.csv --output_suffix pf

Please see the unsupervised clustering ``knitr`` notebook for detailed analysis steps using pre-provided
``heteromotility tools``.
