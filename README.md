pyoptobjects
=============

pyoptobjects is a Python module for optimal object subselection.
The object selection is usefull in computer-based method 
for selecting a subset of similar or diverse object from a large dataset.

Different algorithms have been published and this package aim to collect 
the most important and the most useful algorithms.

This package is distributed under the 3-Clausule BSD license 
and it is currently mantained by Giuseppe Marco Randazzo. 

Voluntary contributions are welcome. :-)

Dependencies
============

The required dependencies to use pyoptobjects is numpy.


Install
=======

To install for all users on Unix/Linux::

  python setup.py build
  sudo python setup.py install

or simply

  easy_install .

Development
===========

GIT
~~~

You can check the latest sources with the command::

  git clone https://github.com/gmrandazzo/pyoptobjects.git
  
~~~

Contributing
~~~~~~~~~~~~~

To contribute you can fork the project, or if you have already forked the project
update to the latest version of pyoptobjects, make the changes and open a Pull Request.

However some recomendation before open a Pull Request:
  * Be sure that your code it's working.
  * Use pylint to check your code. The Global Evaluation rate must be >= 9.0
  * Comment your code with Parameters, Attribute, Return, Notes and References.
  * An example is necessary.
  
Probabily your code will be integrated but some quality and goals have to keep in mind.
