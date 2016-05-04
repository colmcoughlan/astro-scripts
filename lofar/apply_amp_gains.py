from os import system

for i in range(192004,192025+1,3):
	print 'Applying gain amp calibration to L'+str(i)+'_SB182_uv.dppp.MS'
	system('calibrate-stand-alone --parmdb 3C147_gain_amplitude_solutions L'+str(i)+'_SB182_uv.dppp.MS correct_only.bbm > apply_amp_gains_'+str(i)+'.log &')
