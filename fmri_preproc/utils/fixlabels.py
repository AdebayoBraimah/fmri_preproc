# fixlabels.py - Functions for loading/saving FIX/ICA-AROMA label files.
#
# Author: Paul McCarthy <pauldmccarthy@gmail.com>
#
"""This module contains functions for loading/saving FIX/ICA-AROMA label files.

.. autosummary::
    :nosignatures:

    loadLabelFile
    saveLabelFile
    isNoisyComponent
    InvalidLabelFileError
"""


import os.path as op
from typing import List, Optional


def loadLabelFile(
    filename: str,
    includeLabel: List[str] = None,
    excludeLabel: List[str] = None,
    returnIndices: bool = False,
):
    """Loads component labels from the specified file. The file is assuemd
    to be of the format generated by FIX, Melview or ICA-AROMA; such a file
    should have a structure resembling the following::


        filtered_func_data.ica
        1, Signal, False
        2, Unclassified Noise, True
        3, Unknown, False
        4, Signal, False
        5, Unclassified Noise, True
        6, Unclassified Noise, True
        7, Unclassified Noise, True
        8, Signal, False
        [2, 5, 6, 7]


    NOTE:     
        This function will also parse files which only contain a
        component list, e.g.::

            [2, 5, 6, 7]

        The square brackets may or may not be present, i.e. the
        following format is also accepted (this format is generated
        by ICA-AROMA)::

            2, 5, 6, 7

        In this case, the returned melodic directory path will be
        ``None``.  The ``includeLabel`` and ``excludeLabel`` arguments
        allow you to control the labels assigned to included/excluded
        components.


    The first line of the file contains the name of the melodic directory.
    Then, one line is present for each component, containing the following,
    separated by commas:

      - The component index (starting from 1).

      - One or more labels for the component (multiple labels must be
        comma-separated).

      - ``'True'`` if the component has been classified as *bad*,
        ``'False'`` otherwise. This field is optional - if the last
        comma-separated token on a line is not equal (case-insensitive)
        to ``True`` or ``False``, it is interpreted as a component label.

    The last line of the file contains the index (starting from 1) of all
    *bad* components, i.e. those components which are not classified as
    signal or unknown.

    Args:
        filename: Name of the label file to load.
        includeLabel: If the file contains a single line containing a list component indices, this label will be used for the components in the list. Defaults to 'Unclassified noise' for FIX-like files, and 'Movement' for ICA-AROMA-like files.
        excludeLabel: If the file contains a single line containing component indices, this label will be used for the components that are not in the list.  Defaults to 'Signal' for FIX-like files, and 'Unknown' for ICA-AROMA-like files.
        returnIndices: Defaults to ``False``. If ``True``, a list containing the noisy component numbers that were listed in the file is returned.

    Returns:
        Tuple:
            - The path to the melodic directory as specified in the label file
            - A list of lists, one list per component, with each list containing the labels for the corresponding component.
            - If ``returnIndices is True``, a list of the noisy component indices (starting from 1) that were specified in the file.
    """

    signalLabels = None
    filename = op.abspath(filename)

    with open(filename, "rt") as f:
        lines = f.readlines()

    if len(lines) < 1:
        raise InvalidLabelFileError(
            "Invalid FIX classification " "file - not enough lines"
        )

    lines = [l.strip() for l in lines]
    lines = [l for l in lines if l != ""]

    # If the file contains a single
    # line, we assume that it is just
    # a comma-separated list of noise
    # components.
    if len(lines) == 1:

        line = lines[0]

        # if the list is contained in
        # square brackets, we assume
        # that it is a FIX output file,
        # where included components have
        # been classified as noise, and
        # excluded components as signal.
        #
        # Otherwise we assume that it
        # is an AROMA file, where
        # included components have
        # been classified as being due
        # to motion, and excluded
        # components unclassified.
        if includeLabel is None:
            if line[0] == "[":
                includeLabel = "Unclassified noise"
            else:
                includeLabel = "Movement"

        if excludeLabel is None:
            if line[0] == "[":
                excludeLabel = "Signal"
            else:
                excludeLabel = "Unknown"
        else:
            signalLabels = [excludeLabel]

        # Remove any leading/trailing
        # whitespace or brackets.
        line = lines[0].strip(" []")

        melDir = None
        noisyComps = [int(i) for i in line.split(",")]
        allLabels = []

        for i in range(max(noisyComps)):
            if (i + 1) in noisyComps:
                allLabels.append([includeLabel])
            else:
                allLabels.append([excludeLabel])

    # Otherwise, we assume that
    # it is a full label file.
    else:

        melDir = lines[0]
        noisyComps = lines[-1].strip(" []").split(",")
        noisyComps = [c for c in noisyComps if c != ""]
        noisyComps = [int(c) for c in noisyComps]

        # The melodic directory path should
        # either be an absolute path, or
        # be specified relative to the location
        # of the label file.
        if not op.isabs(melDir):
            melDir = op.join(op.dirname(filename), melDir)

        # Parse the labels for every component
        allLabels = []
        for i, compLine in enumerate(lines[1:-1]):

            tokens = compLine.split(",")
            tokens = [t.strip() for t in tokens]

            if len(tokens) < 3:
                raise InvalidLabelFileError(
                    "Invalid FIX classification file - "
                    "line {}: {}".format(i + 1, compLine)
                )

            try:
                compIdx = int(tokens[0])

            except ValueError:
                raise InvalidLabelFileError(
                    "Invalid FIX classification file - "
                    "line {}: {}".format(i + 1, compLine)
                )

            if tokens[-1].lower() in ("true", "false"):
                compLabels = tokens[1:-1]
            else:
                compLabels = tokens[1:]

            if compIdx != i + 1:
                raise InvalidLabelFileError(
                    "Invalid FIX classification file - wrong component "
                    "number at line {}: {}".format(i + 1, compLine)
                )

            allLabels.append(compLabels)

    # Validate the labels against
    # the noisy list - all components
    # in the noisy list should not
    # have 'signal' or 'unknown' labels
    for i, labels in enumerate(allLabels):

        comp = i + 1
        noise = isNoisyComponent(labels, signalLabels)

        if noise and (comp not in noisyComps):
            raise InvalidLabelFileError(
                "Noisy component {} has invalid " "labels: {}".format(comp, labels)
            )

    for comp in noisyComps:

        i = comp - 1
        labels = allLabels[i]
        noise = isNoisyComponent(labels, signalLabels)

        if not noise:
            raise InvalidLabelFileError(
                "Noisy component {} is missing " "a noise label".format(comp)
            )

    if returnIndices:
        return melDir, allLabels, noisyComps
    else:
        return melDir, allLabels


def saveLabelFile(
    allLabels: List[List[str]],
    filename: str,
    dirname: Optional[str] = None,
    listBad: bool = True,
    signalLabels: Optional[List[int]] = None,
):
    """Saves the given classification labels to the specified file. 
    
    The classifications are saved in the format described in the
    :func:`loadLabelFile` method.

    Args:
        allLabels: A list of lists, one list for each component, where each list contains the labels for the corresponding component.
        filename: Name of the file to which the labels should be saved.
        dirname: If provided, is output as the first line of the file. Intended to be a relative path to the MELODIC analysis directory with which this label file is associated. If not provided, a ``'.'`` is output as the first line.
        listBad: If ``True`` (the default), the last line of the file will contain a comma separated list of components which are deemed 'noisy' (see :func:`isNoisyComponent`).
        signalLabels: Labels which should be deemed 'signal' - see the :func:`isNoisyComponent` function.
    """

    lines = []
    noisyComps = []

    # The first line - the melodic directory name
    if dirname is None:
        dirname = "."

    lines.append(dirname)

    # A line for each component
    for i, labels in enumerate(allLabels):

        comp = i + 1
        noise = isNoisyComponent(labels, signalLabels)

        # Make sure there are no
        # commas in any label names
        labels = [l.replace(",", "_") for l in labels]
        tokens = [str(comp)] + labels + [str(noise)]

        lines.append(", ".join(tokens))

        if noise:
            noisyComps.append(comp)

    # A line listing the bad components
    if listBad:
        lines.append("[" + ", ".join([str(c) for c in noisyComps]) + "]")

    with open(filename, "wt") as f:
        f.write("\n".join(lines) + "\n")


def isNoisyComponent(labels: List[str], signalLabels: List[str] = None):
    """Given a set of component labels, returns ``True`` if the component is ultimately classified as noise, ``False`` otherwise.

    Args:
        signalLabels: Labels which are deemed signal. If a component has no labels in this list, it is deemed noise. Defaults to ``['Signal', 'Unknown']``.
    """
    if signalLabels is None:
        signalLabels = ["signal", "unknown"]

    signalLabels = [l.lower() for l in signalLabels]
    labels = [l.lower() for l in labels]
    noise = not any([sl in labels for sl in signalLabels])

    return noise


class InvalidLabelFileError(Exception):
    """Exception raised by the :func:`loadLabelFile` function when an attempt is made to load an invalid label file.
    """

    pass
