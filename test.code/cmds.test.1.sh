
func=$(realpath data/sub-198/func/sub-198_task-rest_dir-PA_run-01_bold.nii.gz)

ap=$(realpath data/sub-198/func/sub-198_task-rest_dir-AP_run-01_sbref.nii.gz)   # Reversed phase-encoded b0
pa=$(realpath data/sub-198/func/sub-198_task-rest_dir-PA_run-01_sbref.nii.gz)   # Phase-encoded b0

echo "0 1 0 0.05" >> param.acqp
echo "0 -1 0 0.05" >> param.acqp

fslmerge -t b0s ${pa} ${ap}

z_dim=$(fslval ${ap} dim3)
nvols=$(fslval ${func} dim4)

for (( i = 0; i < ${nvols}; i++ )); do
  echo 1 >> frame.idx
done

# check if even
#   if even - do nothing
#   if odd - remove bottom slice

fslroi b0s.nii.gz b0s_even 0 -1 0 -1 1 $((z_dim - 1))

topup \
--imain=b0s_even \
--datain=param.acqp \
--config=${FSLDIR}/etc/flirtsch/b02b0.cnf \
--out=suscept_corr_B0 \
--fout=suscept_field_Hz \
--iout=unwarped_B0s \
--scale=1 --verbose

fslsplit b0s_even b0.tmp_ -t

mv b0.tmp_*0.nii.gz B0_PA.nii.gz
mv b0.tmp_*1.nii.gz B0_AP.nii.gz

applytopup \
--imain=B0_PA.nii.gz,B0_AP.nii.gz \
--datain=param.acqp \
--inindex=1,2 \
--topup=suscept_corr_B0 \
--method=lsr \
--out=hifi \
--verbose

bet hifi.nii.gz hifi_brain -m -R -f 0.25

# generate bvals & bvecs
#   * See jupyter notebook

# Eddy current/distortion/s2v correction
#   NOTE: Test s2v parameters on cluster

module load fsl/6.0.4 
module load cuda/9.1

fslroi ${func} func_even 0 -1 0 -1 1 $((z_dim - 1))

# bsub -M 10000 -W 1000 -n 1 -J "eddy_test" eddy \
bsub -M 10000 -W 1000 -n 1 -J "eddy_test" -R "rusage[gpu=1]" -q gpu-v100 eddy_cuda \
--imain=func_even \
--mask=hifi_brain_mask.nii.gz \
--index=file.idx \
--acqp=file_functional.acqp \
--bvecs=file.bvec \
--bvals=file.bval \
--out=eddy_corr_fmri \
--topup=suscept_corr_B0 \
--estimate_move_by_susceptibility \
--data_is_shelled \
--very_verbose \
--b0_only \
--dont_mask_output \
--s2v_interp=10 \
--s2v_fwhm=0 \
--slspec=file.slice.order \
--mporder=10 \
--s2v_interp=trilinear \
--s2v_lambda=1
# --residuals       # Not relevant for B0s
# --repol           # Not relevant for B0s
# --cnr_maps        # Not relevant for B0s
