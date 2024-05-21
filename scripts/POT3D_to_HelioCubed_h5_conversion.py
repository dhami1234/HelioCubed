import h5py
import numpy as np
from datetime import datetime as dt
import time
from glob import glob

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

#Initialize
r0 = 0.1
num_components = 8
time = toYearFraction(dt(2020, 2, 2, 22, 0, 0))
phys_domain = [0.1, 0, 0, 0.1, 6.28319, 3.14159]
# Out_file = "POT3D360x180_2020020222_fixeddtheta_Ron.hdf5"
Out_file = "test.h5"
step_const = [0, 1, 0]

address_of_files = '/Users/talwinder/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/POT3D_to_HCubed_h5/Time_independent/'
# address_of_files = '/Users/talwinder/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/POT3D_to_HCubed_h5/Time_dependent/'

files = glob(address_of_files+'rho_upper_corona_Np*.h5')
files=sorted(files)
rho_filename = address_of_files+"rho_upper_corona_Np.h5"
V_filename = address_of_files+"vr_upper_corona_kms.h5"
T_filename = address_of_files+"t_upper_corona_K.h5"
B_filename = address_of_files+"br_upper_corona_Gauss.h5"

#If n time dependent BCs availaible, they should be numbered like rho_upper_corona_Np1, ..., rho_upper_corona_Npn

for i in range(len(files)):
    print("Converting file:",i+1)
    if (len(files) > 1):
        rho_filename = address_of_files+"rho_upper_corona_Np"+str(i+1)+".h5"
        V_filename = address_of_files+"vr_upper_corona_kms"+str(i+1)+".h5"
        T_filename = address_of_files+"t_upper_corona_K"+str(i+1)+".h5"
        B_filename = address_of_files+"br_upper_corona_Gauss"+str(i+1)+".h5"
    with h5py.File(rho_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_rho = list(f[a_group_key])
        data_rho = np.array(data_rho)
        dim_siz1 = data_rho.shape[0]
        dim_siz2 = data_rho.shape[1]
        data_rho = data_rho.flatten('F')
        b_group_key = list(f.keys())[1]
        theta = list(f[b_group_key])
        theta = np.array(theta)
        theta = theta.flatten('F')
        c_group_key = list(f.keys())[2]
        phi = list(f[c_group_key])
        phi = np.array(phi)
        phi = phi.flatten('F')

    with h5py.File(V_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        print(list(f.keys()))
        data_Vr = list(f[a_group_key])
        data_Vr = np.array(data_Vr)
        data_Vr = data_Vr.flatten('F')
        data_Vr = data_Vr*1e5  #make it cm/s
        data_Vp = data_Vr*0.0
        data_Vt = data_Vr*0.0

    with h5py.File(T_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_T = list(f[a_group_key])
        data_T = np.array(data_T)
        data_T = data_T.flatten('F')   
        data_T = data_T*1e-6  #make it Mega Kelvin

    data_P = 2.0*data_rho*1.3806505e-16*data_T/1e-12   # in pico Dyne

    with h5py.File(B_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_Br = list(f[a_group_key])
        data_Br = np.array(data_Br)
        data_Br = data_Br.flatten('F')
        data_Br = data_Br*1e6  #make it micro Gauss
        data_Bp = data_Br*0.0
        data_Bt = data_Br*0.0 

    data_final = np.concatenate((data_rho,data_Vr,data_Vp, data_Vt, data_P, data_Br, data_Bp, data_Bt))
    dtheta = np.full(dim_siz2, np.pi/dim_siz2)
    if (i==0):
        with h5py.File(Out_file, "w") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final,dtype='float64')
            data_file["data"+str(i)].attrs["time"] = i
            data_file.attrs['domain'] = [dim_siz1,dim_siz2]
            data_file.attrs['num_components'] = num_components
            data_file.attrs['num_datasets'] = len(files)
            data_file.attrs['time'] = time
            data_file.attrs['r0'] = r0
            data_file.attrs['datasets_time'] = range(len(files))
            grp = data_file.create_group("geometry")
            grp.attrs['phys_domain'] = phys_domain
            grp.attrs['step_const'] = step_const
            grp.create_dataset("dtheta", data=dtheta,dtype='float64')
            grp.create_dataset("theta", data=theta,dtype='float64')
            grp.create_dataset("phi", data=phi,dtype='float64')
    if (i!=0):
        with h5py.File(Out_file, "a") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final,dtype='float64')
            data_file["data"+str(i)].attrs["time"] = i