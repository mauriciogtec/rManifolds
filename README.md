rManifolds 
================


SnapPy/Python module for creating random manifolds from knots and braids. Python package index: https://pypi.python.org/pypi/rManifolds.

Currently, there are only three variants supported:
* *Gaussian knots*: uses independent normal distribution to simulate the deplacement in a stochastic processes than begins and ends in the same place.
* *Random uniform knots*: generates a stochastic processes in which the position of the next vertex is selected uniformly from a fixed confined space.
* *Random braids*: uses a discrete probability distribution in the generators of the braid group to create a random braid of a fixed number of strands. Currently, it always uses a uniform distribution on the generators, although general distributions will be added.


We know describe each of these methods.

```python
from rManifolds import *
```

<h2> Gaussian knots </h2>

A gaussian random walk of variance <b>sigma^2</b> is defined as a random walk in which the deplacements in each direction are independent and normally distribute with mean 0 and variance <b>sigma^2</b>. Given a fixed number of steps, this stochastic processes is conditioned so that the last point is the same as the departing point, obtaining a knot. The random walk is modeled in R^3 and then scaled to fit nicely into the SnapPy's standard grid of 500x500 inside the PLink editor.

This model of a random knot has found some applications for modeling protein molecules [1]. As the number of vertices goes to infinity, these knots usually turn out to be satellite knots due to the formation of clusters [2].

```python
In [1]: n = 1000 
In [2]: obj = rGaussianKnot(n, filename='example1', sigma=1) 
A file 'example1.lnk' has been created in the current working directory.
```
The default value of sigma is 1 when it is not specified. If a filename is not specified then the .lnk object is not created. The filename could is allowed to be a simple string name, in which case the file is saved in the current working directory, or it could be a long directory address, in which case it is saved in the specified address. Once the object has been created, we can recover the manifold and the filename from the object. Methods for manifold classes are available.

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

The .lnk file created can be opened using the SnapPy interface and the Plink editor. The following image is an example of Gaussian knot with a 1000 vertices.

<img src="https://github.com/mauriciogtec/rManifolds/blob/master/screenshots/GaussianKnot.png?raw=true" alt="alt text" width="250" height="250">

<h2> Random uniform knots </h2>
This is the simplest model of rando knot. Random points in a cube are generated and and joint to form a knot (the initial and the last point are joint together). A simple system of equation can be solved to create a planar diagram of the knot. Both the filename and the manifold can be recovered and manipulated as in the case of Gaussian knots.

```python
In [1]: M = rUnifKnot(100)
```

Here is a screenshot of resulting knot:

<img src="https://github.com/mauriciogtec/rManifolds/blob/master/screenshots/rUnifKnot.png?raw=true" alt="alt text" width="250" height="250">

<h2> Random braids </h2>

A different way of generating random manifolds is by specifying a probability distribution in the Braid group of n strands. In the current version a uniform distribution is given over the generators of the braid group. To create a Braid object one needs two parameters, one for the numbers of strands, and another one for the number of crossings in the link diagram. The result is sometimes a knot and sometimes a link. Recently, Ma [3] proved that as the length of the walk goes to infinity, the probability that the resulting link is hyperbolic goes to one.

From the created random braid one can recover the filename (when specified) and the manifold object as in the above examples (see Gaussian knots). For random braids, one can also ask for the number of link components.

```python
In [1]: M = rBraid(n=5,k=40)
In [2]: M.num_components
Out[2]: 3
In [3]: M.manifold.homology()
Out[3]: Z + Z + Z
In [4]: M.manifold.volume()
Out[4]: 26.635490252
```
Here is a screenshot of a random braid:

<img src="https://github.com/mauriciogtec/rManifolds/blob/master/screenshots/rBraid.png?raw=true" alt="alt text" width="250" height="250">

<h3> References </h3>
[1] Diao et al. (1994). "On random knots". &lt;em&gt; Random knotting and linking &lt;/em&gt;. World Scientific Publishing Co.

[2] Adams, C. (2004). "Hyperbolic knots". &lt;em&gt; Handbook of Knot theory &lt;/em&gt;. http://arxiv.org/pdf/math/0309466.pdf

[3] Ma, J. (2014). "The closure of a random braid is a hyperbolic link". &lt;em&gt; Proceedings of the AMS &lt;/em&gt; 142. pp. 695-701.
