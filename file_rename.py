import glob
import os

current_basename = 'test'
new_basename     = 'guider'

files = glob.glob(current_basename+'*.fits')
for f in files: 
	new_name = f.replace(current_basename,new_basename)
	os.rename(f,new_name)
