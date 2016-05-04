export o_dir='/mnt/home_cr/coughlan/lofar/data/data/TTAU_SB360_369_sagecal'
export working_dir='.'
export script_dir='/mnt/home_cr/coughlan/lofar/scripts'
export subband='360'
export num_subbands='10'

export l_start='192737'
export l_end='192758'

export nthreads='4'
export gsmskymodel='/mnt/home_cr/coughlan/lofar/skymodels/GSMSkymodel'
#export parmpath='/mnt/home_cr/coughlan/lofar/parmexportcal_plus_clock/parmexportcal_plus_clock'
export parmpath='/home/coughlan/phi_stuff/parmexportcal_plus_clock/parmexportcal_plus_clock'

echo "Starting initial amp and phase cal..."

python ../scripts/average.py

echo "Initial amp and phase cal complete."
