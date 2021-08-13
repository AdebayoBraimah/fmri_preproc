# -*- coding: utf-8 -*-
"""fMRI preprocessing pipeline - main pipeline.
"""
# Pipeline overview
# 
# 1. Prepare fieldmap
#   a. Topup [mc]
# 2. Motion correction, distortion correction
#   a. eddy [mc]
#   b. mcflirt [mc]
# 3. Registration (linear, and non-linear)
#   a. flirt ( + mcflirt) [reg]
#       i. sym-link field map files from [mc]
# 4. ICA
#   a. meldic
# 5. Denoise
#   a. FIX
# 6. QC
#   a. PENDING

class Pipeline:
    def __init__(self, 
                 outdir: str
                ) -> None:
        """Constructor for pipeline class.
        """
        pass

    def prepare_fieldmap(self) -> None:
        pass

    def mcdc(self) -> None:
        pass

    def standard(self) -> None:
        pass

    def ica(self) -> None:
        pass

    def denoise(self) -> None:
        pass

    def report(self) -> None:
        pass

    def qc(self) -> None:
        pass

    def pre_mcdc(self) -> None:
        """Perform all pre-motion-and-distortion-correction stages of the preprocessing pipeline."""
        pass

    def post_mcdc(self) -> None:
        """Perform all post-motion-and-distortion-correction stages of the preprocessing pipeline."""
        pass

    def run_all(self) -> None:
        """Perform all stages of the preprocessing pipeline."""
        pass
