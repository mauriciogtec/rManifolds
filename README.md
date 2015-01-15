rManifolds v0.1
================


SnapPy/Python module for creating random manifolds from knots and braids. Currently, there are only three variants supported:
* *Gaussian knots*: uses independent normal distribution to simulate the deplacement in a stochastic processes than begins and ends in the same place.
* *Random uniform knots*: generates a stochastic processes in which the position of the next vertex is selected uniformly from a fixed confined space.
* *Random braids*: uses a discrete probability distribution in the generators of the braid group to create a random braid of a fixed number of strands. Currently, it always uses a uniform distribution on the generators, although general distributions will be added.

We know describe each of these methods.

```python
from rManifolds import *
```

<h2> Gaussian knots </h2>

A gaussian random walk of variance <b>sigma^2</b> is defined as a random walk in which the deplacements in each direction are independent and normally distribute with mean 0 and variance <b>sigma^2</b>. Given a fixed number of steps, this stochastic processes is conditioned so that the last point is the same as the departing point, obtaining a knot. The random walk is modeled in R^3 and then scaled to fit nicely into the SnapPy's standard grid of 500x500 inside the PLink editor.

This model of a random knot has found some applications for modeling protein molecules [1]. As the number of vertices goes to infinity, these knots usually turn out to be satellite knots due to the formation of clusters [3].

```python
In [1]: n = 1000 # a parameter for the number of vertices in the knot.
In [2]: obj = rGaussianKnot(n, filename='example1', sigma=1) # this will store as a rManifold object.
A file named "example1.lnk" was created in the current working directory.
```
The default value of sigma is 1. If a filename is not specified then the .lnk object is not created. The filename could be just a string, in which case the file is saved in the current working directory, or it could be a long address, in which case it is saved in the specified address. Once the object has been created, we can recover the manifold and the filename from the object. Methods for the objects of type manifolds are available.

```python
In [3]: obj.filename
Out[3]: '/home/mauriciogtec/example1.lnk'
In [4]: obj.manifold
Out[4]: unnamed link(0,0)
In [5]: obj.manifold.homology()
Out[5]: Z
In [6]: obj.manifold.volume()
Out[6]: 1.0 E-05
```
