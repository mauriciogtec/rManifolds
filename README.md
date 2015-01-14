rManifolds
==========

SnapPy/Python module for creating random manifolds.

```python
from rManifolds import *
```

<h2> Gaussian knots </h2>

A gaussian random walk of variance <b>sigma^2</b> is defined as a random walk in which the deplacements in each direction are independent and normally distribute with mean 0 and variance <b>sigma^2</b>. Given a fixed number of steps, this stochastic processes is conditioned so that the last point is the same as the departing point, obtaining a knot. The random walk is modeled in R^3 and then scaled to fit nicely into the SnapPy's standard grid of 500x500 inside the PLink editor.

This model of a random knot has found some applications for modeling protein molecules [1]. As the number of vertices goes to infinity, these knots usually turn out to be satellite knots due to the formation of clusters [3].

