from pyrap import tables as pt
import numpy as np

fn='L192737_SB350_uv.dppp.flagged.shifted.super.scaled.MS'
station_list = [0,1,2,3,4,5,20,21]
data_col = 'DATA'
int_int_factor=100000.0
int_crs_factor = np.sqrt(int_int_factor)

print('Rescaling visibilities of stations '+str(station_list)+' in '+data_col+' column of '+fn)

t = pt.table(fn , readonly=False)

t1 = t.query('(ANTENNA1 in '+str(station_list)+')&&(ANTENNA2 in '+str(station_list)+')')
data = t1.getcol(data_col)
print('Selecting INT-INT baselines only')
print('Detected '+str(len(data))+' rows')
if(len(data)>0):
	print('Scaling selected station fluxes by a factor of '+str(int_int_factor))
	data = np.multiply(data , int_int_factor)
	t1.putcol(data_col,data)
	t.flush()
else:
	print('No stations detected.')
	print('Error! - this script assumes interational stations')

t1 = t.query('(ANTENNA1 in '+str(station_list)+')&&!(ANTENNA2 in '+str(station_list)+')')
data = t1.getcol(data_col)
print('Selecting INT-[TR]S baselines only')
print('Detected '+str(len(data))+' rows')
if(len(data)>0):
	print('Scaling selected station fluxes by a factor of '+str(int_crs_factor))
	data = np.multiply(data , int_crs_factor)
	t1.putcol(data_col,data)
	t.flush()
else:
	print('No stations detected.')

t1 = t.query('!(ANTENNA1 in '+str(station_list)+')&&(ANTENNA2 in '+str(station_list)+')')
data = t1.getcol(data_col)
print('Selecting [TR]S-INT baselines only')
print('Detected '+str(len(data))+' rows')
if(len(data)>0):
	print('Scaling selected station fluxes by a factor of '+str(int_crs_factor))
	data = np.multiply(data , int_crs_factor)
	t1.putcol(data_col,data)
	t.flush()
else:
	print('No stations detected.')

