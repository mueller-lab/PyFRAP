#Import setuptools
from setuptools import setup

#We need those two for doing some file permission magic later
import os
import platform

#Overwrite setuptools install to be able to set file permissions
from setuptools.command.install import install
from distutils import log 
from setuptools.command.install_scripts import install_scripts

#Override of install that allows us to set file permissions for meshfiles and configurations
#Idea taken from http://stackoverflow.com/questions/5932804/set-file-permission-in-setup-py-file (thanks a bunch!)
class OverrideInstall(install):

	def run(self):
		
		#Grab user ID and group ID of actual user
		uid=pwd.getpwnam(os.getlogin())[2]
		gid=pwd.getpwnam(os.getlogin())[3]
		
		#Mode for files (everyone can read/write/execute. This is somewhat an overkill, but 0666 seems somehow not to work.)
		mode=0777
		
		#Run setuptools install
		install.run(self) 
		
		#Print log info
		log.info("Overriding setuptools mode of scripts ...")
		
		#Overwrite file permissions
		for filepath in self.get_outputs():
			
			if "meshfiles" in filepath or "configurations" in filepath:
				
				#Change permissions for folder containing files
				folderpath=os.path.dirname(os.path.realpath(filepath))
				self.changePermissions(folderpath,uid,gid,mode)
				
				#Change permissions of file
				self.changePermissions(filepath,uid,gid,mode)
				
				#Make some more necessary data folders
				if folderpath.endswith("meshfiles"):	
					self.makeAdditionalDataFolders(folderpath,"field",uid,gid,mode)
					self.makeAdditionalDataFolders(folderpath,"field/custom",uid,gid,mode)
	
	def changePermissions(self,filepath,uid,gid,mode):
		ret=True
		try:
			os.chown(filepath, uid, gid)
			log.info("Changing ownership of %s to uid:%s gid %s" %(filepath, uid, gid))
		except:
			log.info("Was not able to change ownership of file %s" %(filepath))
			ret=False

		try:
			log.info("Changing permissions of %s to %s" %(filepath, oct(mode)))
			os.chmod(filepath, mode)
		except:	
			log.info("Was not able to change file permissions of file %s" %(filepath))
			ret=False
		return ret
		
	def makeAdditionalDataFolders(self,folder,fn,uid,gid,mode):
		if not folder.endswith("/"):
			folder=folder+"/"
			
		if os.path.isdir(folder+fn):
			return False
		else:
			try:
				os.mkdir(folder+fn)
				self.changePermissions(folder+fn,uid,gid,mode)
				return True
			except:
				log.info("Unable to create folder %s" %(folder+fn))
				return False
				
#Define setup

#Check which operating system: If Windows, do not override installation procedure
if platform.system() in ["Windows"]:

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
	)

#If OS is not Windows, we will overwrite the install procedure to set the file permissions for data files
else:
	#Import this only here, since pwd only exists for Unixoids.
	import pwd
	
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
