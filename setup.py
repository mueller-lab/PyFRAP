#Import setuptools
from setuptools import setup

#We need those two for doing some file permission magic later
import os
import platform
import shutil

#Overwrite setuptools install to be able to set file permissions
from setuptools.command.install import install
from distutils import log 
from setuptools.command.install_scripts import install_scripts

#Override of install that allows us to set file permissions for meshfiles and configurations
#Idea taken from http://stackoverflow.com/questions/5932804/set-file-permission-in-setup-py-file (thanks a bunch!)
class OverrideInstall(install):

	def run(self):
		
		#Try to download gmsh
		#self.downloadGmsh()
				
		#Run setuptools install
		install.run(self) 
		
		#Print log info
		log.info("Overriding setuptools mode of scripts ...")
		
		#Add Data and edit file permissions
		self.addData()
			
	def addData(self):
		
		if platform.system() not in ["Windows"]:
			import pwd
			
			#Grab user ID and group ID of actual user
			uid=pwd.getpwnam(os.getlogin())[2]
			gid=pwd.getpwnam(os.getlogin())[3]
			
			#Mode for files (everyone can read/write/execute. This is somewhat an overkill, but 0666 seems somehow not to work.)
			mode=0777
		
		#Overwrite file permissions
		for filepath in self.get_outputs():
			
			if platform.system() not in ["Windows"]:
			
				if "meshfiles" in filepath or "configurations" in filepath or "executables" in filepath:
					
					#Change permissions for folder containing files
					folderpath=os.path.dirname(os.path.realpath(filepath))
					self.changePermissions(folderpath,uid,gid,mode)
					
					#Change permissions of file
					self.changePermissions(filepath,uid,gid,mode)
					
					#Make some more necessary data folders
					if folderpath.endswith("meshfiles"):
						self.makeAdditionalDataFolders(folderpath,"field",uid,gid,mode)
						self.makeAdditionalDataFolders(folderpath,"field/custom",uid,gid,mode)
						
					if folderpath.endswith("configurations"):
						self.makeAdditionalDataFolders(folderpath,"macros",uid,gid,mode)
			
			#Add gmsh into paths.default if download was successful
			#if 'paths.default' in filepath:
				#if self.gmshDownloaded:
					#self.setGmshPath(filepath)
					#if platform.system() not in ["Windows"]:
						#self.changePermissions(filepath,uid,gid,mode)
				
	def downloadGmsh(self):
		
		#Define gmshVersion (might need to update this line once in a while)
		gmshVersion='2.12.0'
		
		#Flag to see if gmsh DL went through
		self.gmshDownloaded=False
		
		#Make executables folder if it doesn't exist yet
		try:
			os.mkdir('pyfrp/executables')	
		except OSError:
			log.info('Was not able to create directore pyfrp/executables')
		
		#Get fileList before
		filesBefore=os.listdir('.')
		
		#Try to import wget
		try: 
			import wget
			
			#Get Architecture
			arch=platform.architecture()[0].replace('bit','')
			
			#For windows
			if platform.system() in ["Windows"]:
				
				#Download Gmsh
				url='http://gmsh.info/bin/Windows/gmsh-'+gmshVersion+'-Windows'+arch+'.zip'
				folderFn=wget.download(url)
				fnDL=str(folderFn)
				print
				
				#Decompress
				import zipfile 
				with zipfile.ZipFile(folderFn) as zf:
					zf.extractall()
					
				folderFn='gmsh-'+gmshVersion+'-Windows'	
				
				self.gmshPath='executables/gmsh/bin/gmsh.exe'
				
			elif platform.system() in ["Linux"]:
				
				#Download Gmsh
				url='http://gmsh.info/bin/Linux/gmsh-'+gmshVersion+'-Linux'+arch+'.tgz'
				
				folderFn=wget.download(url)
				fnDL=str(folderFn)
				print
				
				#Decompress
				import tarfile
				with tarfile.open(folderFn,mode='r:gz') as zf:
					zf.extractall()
				
				folderFn='gmsh-'+gmshVersion+'-Linux'
				
				self.gmshPath='executables/gmsh/bin/./gmsh'
				
			elif platform.system() 	in ["Darwin"]:
				
				#Download Gmsh
				url='http://gmsh.info/bin/MacOSX/gmsh-'+gmshVersion+'-MacOSX'+'.dmg'
				folderFn=wget.download(url)
				
				fnDL=str(folderFn)
				
				#Mount dmg file (Here the user need to read through LICENSE, don't know how to fix this)
				os.system('hdiutil attach '+folderFn)
				folderFn=folderFn.replace('.dmg','')
				try:
					os.mkdir(folderFn)
				except OSError:
					pass
				
				cwd=os.getcwd()
				#Copy gmsh executable to cwd
				os.system('cp -rv /Volumes/'+folderFn+'/Gmsh.app/Contents/MacOS/bin/ '+ cwd)
				os.system('cp -rv /Volumes/'+folderFn+'/Gmsh.app/Contents/MacOS/share/ '+ cwd)
				
				#Unmount gmsh
				os.system('hdiutil detach /Volumes/'+folderFn+'/')
		
				self.gmshPath='executables/gmsh/bin/./gmsh'
					
			#Copy file to pyfrp/executables/
			try:
				shutil.rmtree('pyfrp/executables/gmsh/')
			except:
				pass
			shutil.copytree(folderFn+"/",'pyfrp/executables/gmsh')
			
			#Remove downloaded files
			os.remove(fnDL)
			
			#Get fileList before
			filesAfter=os.listdir('.')
			
			#Difference between files
			filesDiff=list(set(filesAfter)-set(filesBefore))
			try:
				shutil.rmtree(filesDiff[0])
			except:
				pass
			
			#Remove files
			log.info("Installed gmsh to "+ self.gmshPath)
			
			#Set Flag=True
			self.gmshDownloaded=True
			
		except ImportError:
			log.info("Cannot find wget, will not be downloading gmsh. You will need to install it later manually")
	
	def setGmshPath(self,fn):
		
		#Make backup of default path file
		shutil.copy(fn,fn+'_backup')
		
		#Get filepath to PyFRAP
		fnPyfrp=fn.replace('configurations/paths.default','')
		
		#Open file and enter new gmsh bin
		with open(fn,'rb') as fPath:
			with open(fn+"_new",'wb') as fPathNew:
				for line in fPath:
					if line.strip().startswith('gmshBin'):
						ident,path=line.split('=')
						path=path.strip()
						lineNew=ident+"="+fnPyfrp+self.gmshPath
						fPathNew.write(lineNew+'\n')
					else:
						fPathNew.write(line)
			
		#Rename file
		shutil.move(fn+'_new',fn)
		
	def changePermissions(self,filepath,uid,gid,mode,fullOutput=False):
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

import os

if os.environ.get('READTHEDOCS', None) == 'True':
	
	print "Installing on RTD, will not overwrite install command."
	
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
		zip_safe=False
		)

else:

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
