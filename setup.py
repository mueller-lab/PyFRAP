

#Import setuptools
from setuptools import setup

#We need those two for doing some file permission magic later
import os
import pwd

#Overwrite setuptools install to be able to set file permissions
from setuptools.command.install import install
from distutils import log 
from setuptools.command.install_scripts import install_scripts

#Override of install that allows us to set file permissions for meshfiles and configurations
class OverrideInstall(install):

	def run(self):
		
		#Grab user ID and group ID of actual user
		uid=pwd.getpwnam(os.getlogin())[2]
		gid=pwd.getpwnam(os.getlogin())[3]
		
		#Mode for files (owner can read/write, everyone else read)
		mode=0644
		
		#Run setuptools install
		install.run(self) 
		
		#Print log info
		log.info("Overriding setuptools mode of scripts ...")
		
		#Overwrite file permissions
		for filepath in self.get_outputs():
	
			if "meshfiles" in filepath or "configurations" in filepath:
							
				try:
					os.chown(filepath, uid, gid)
					log.info("Changing ownership of %s to uid:%s gid %s" %(filepath, uid, gid))
				except:
					log.info("Was not able to change ownership of file %s" %(filepath))
		
				try:
					log.info("Changing permissions of %s to %s" %(filepath, oct(mode)))
					os.chmod(filepath, mode)
				except:	
					log.info("Was not able to change file permissions of file %s" %(filepath))

#Define setup
setup(name='pyfrp',
      version='1.0',
      description='PyFRAP: A Python based FRAP analysis tool box',
      url='https://github.com/alexblaessle/PyFRAP',
      author='Alexander Blaessle',
      author_email='alexander.blaessle@tuebingen.mpg.de',
      license='GNU GPL Version 3',
      packages=['pyfrp','pyfrp.modules','pyfrp.subclasses','pyfrp.gui'],
      package_dir={'pyfrp': 'pyfrp',
		   'pyfrp.modules': 'pyfrp/modules',
		   'pyfrp.gui' : 'pyfrp/gui'
		   },
      #package_data = {'pyfrp':['meshfiles','configurations']},
      include_package_data=True,
      classifiers= [
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Biophysics/FRAP :: Analysis/Visualization',
        'Intended Audience :: Science/Research'
        ],
      platforms=['ALL'],
      keywords=["FRAP", "fluorescence",'recovery','after','photobleaching','reaction-diffusion','fitting'
              ],
      zip_safe=False,
      cmdclass={'install': OverrideInstall} #Need this here to overwrite our install
      )
      


