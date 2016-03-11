import os, site, sys

#Remember old path
old_user_base = site.USER_BASE


#Rename Paths
site.USER_BASE = '/fml/ag-mueller_share/PyFRAP_and_IdealFRAPData/python-packages/'
site.USER_SITE = site.USER_SITE.replace(old_user_base, site.USER_BASE)

#Make site directory
if not os.path.exists(site.USER_SITE):
	os.makedirs(site.USER_SITE)
	site.addsitedir(site.USER_SITE, known_paths=[])
	
#Add to PATH	
for i, path in enumerate(sys.path[:]):
	if old_user_base in path:
		sys.path.remove(path)
		new_path = path.replace(old_user_base, site.USER_BASE)
		if os.path.exists(new_path):
			sys.path.insert(i, new_path)
	
if site.USER_SITE not in sys.path:
	sys.path.insert(0, site.USER_SITE)