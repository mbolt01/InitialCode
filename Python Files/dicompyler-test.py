# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 19:09:31 2017

@author: mb22
"""

from dicompylercore import dicomparser, dvh, dvhcalc
import os

#%%
## import the dicom dose and structures
user = os.getlogin()
root_path = r'C:\Users\\'
folder_path = '\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\DoseCubes'
file_dose = '\PARSPORT_1001\RD.1.3.6.1.4.1.25111.0000.2.000005.20100628.170040.dcm'
file_struct = '\PARSPORT_1001\RS.1.3.6.1.4.1.25111.0000.4.000000.20100628.170040.dcm'
#dicom_dose = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\DoseCubes\PARSPORT_1001\RD.1.3.6.1.4.1.25111.0000.2.000005.20100628.170040.dcm'
#dicom_struct = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\DoseCubes\PARSPORT_1001\RS.1.3.6.1.4.1.25111.0000.4.000000.20100628.170040.dcm'

dicom_dose = root_path +user + folder_path + file_dose
dicom_struct = root_path +user + folder_path + file_struct
print(dicom_dose)
print(dicom_struct)

rtdose = dicomparser.DicomParser(dicom_dose)
#%%
structures = dicomparser.DicomParser(dicom_struct).GetStructures()
for i in range(1,5):
    print(structures[i])
#%%
## calcualte dvh for given structure ID
calcdvh = dvhcalc.get_dvh(dicom_struct, dicom_dose, 28)
## DVH stats
calcdvh.max, calcdvh.min, calcdvh.D2cc, calcdvh.mean
## decription of DVH.
## might be able to edit the code to give more values
## OR look at the function to see if can get the individual values myself directly
calcdvh.describe()

#%%
## thsi might give me ability to get stats....
stat = calcdvh.D90
print(stat)

#%%

dvh.DVH.from_data()