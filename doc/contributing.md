# Contributing to ``fmri_preproc``
-----------------------------------

> This document is adapted from [https://git.fmrib.ox.ac.uk/fsl/fslpy/-/blob/master/doc/contributing.rst](https://git.fmrib.ox.ac.uk/fsl/fslpy/-/blob/master/doc/contributing.rst)        
> Created by Paul McCarthy

**NOTE**: This document is a work in progress.

Development model
-----------------


* All development occurs on separate ``dev`` branches.

* All changes to the master branch occur via merge requests. Individual
  developers are free to choose their own development workflow in their own
  repositories.


Commit messages
---------------


To aid readability, all commit messages should be prefixed with one or more of
the following labels (this convention has been inherited from [nibabel](https://github.com/nipy/nibabel):

  * *BF*  : bug fix
  * *RF*  : refactoring
  * *ENH*:  enhancement/new feature
  * *BW*  : addresses backward-compatibility
  * *OPT* : optimization
  * *BK*  : breaks something and/or tests fail
  * *PL*  : making pylint happier
  * *DOC* : for all kinds of documentation related commits
  * *TEST*: for adding or changing tests
  * *MNT* : for administrative/maintenance changes
  * *CI*  : for continuous-integration changes


Version number
--------------


The ``fmri_preproc`` version number roughly follows [semantic versioning](http://semver.org/) 
rules, so that dependant projects are able to perform
compatibility testing.  The full version number string consists of three
numbers::

    major.minor.patch

- The ``patch`` number is incremented on bugfixes and minor
  (backwards-compatible) changes.

- The ``minor`` number is incremented on feature additions and/or
  backwards-compatible changes.

- The ``major`` number is incremented on major feature additions, and
  backwards-incompatible changes.


The version number in the ``master`` branch should be of the form
``major.minor.patch.dev0``, to indicate that any releases made from this
branch are development releases (although development releases are not part of
the release model).