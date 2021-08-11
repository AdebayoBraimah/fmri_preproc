#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Adebayo Braimah
# Imaging Research Center
# 
# Copyright 2021 Cincinnati Children's Hospital Medical Center
# 
# TODO: <LICENCSE HERE>
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""This :mod: `fmri_preproc` package contains modules for the `fmri_preproc` fMRI pipeline. 
It is split into the following sub-packages:

.. autosummary::

    fsl.func
    fsl.util
"""
import os
name: str = "fmri_preproc"

_version_file: str = os.path.abspath(
    os.path.join('resources','version.txt')
)

with open(_version_file,'r') as f:
    file_contents: str = f.read()
    _version: str = file_contents.strip("\n")
    f.close()

__author__          = "Adebayo Braimah"
__credits__         = ["Adebayo Braimah",
                       "Cincinnati Children's Hospital Medical Center", 
                       "Imaging Research Center", 
                       "CCHMC Dept. of Radiology"]
__license__         = "<LICENSE HERE>"
__version__         = _version
__maintainer__      = "Adebayo Braimah"
__email__           = "adebayo.braimah@cchmc.org"
__status__          = "Development"
