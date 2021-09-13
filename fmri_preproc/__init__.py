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

    func
    util
"""
import os
name: str = "fmri_preproc"


# Define constants
RESOURCEDIR: str = os.path.join(os.path.dirname(__file__), 'resources')
ATLASDIR: str = os.path.join(RESOURCEDIR,'atlases')
HTMLDIR: str = os.path.join(RESOURCEDIR,'html')
FIXDATADIR: str = os.path.join(RESOURCEDIR,'fix')
GROUP_MAP_DIR: str = os.path.join(RESOURCEDIR,'group_maps')
GROUP_QC_DIR: str = os.path.join(RESOURCEDIR,'group_qc')
SETTINGS: str = os.path.join(RESOURCEDIR,'settings.default.json')


_version_file: str = os.path.abspath(
    os.path.join(RESOURCEDIR,'version.txt')
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
