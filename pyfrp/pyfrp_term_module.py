import colorama
import PyQt4.QtGui as QtGui

def printWarning(txt):

	print(colorama.Fore.YELLOW + "WARNING:") + colorama.Fore.RESET + txt

def printError(txt):
	
	print(colorama.Fore.RED + "ERROR:") + colorama.Fore.RESET + txt

def printNote(txt):

	print(colorama.Fore.GREEN + "NOTE:") + colorama.Fore.RESET + txt
	
def printDict(dic):
	for k in dic.keys():
		print k , ' = ' , dic[k]
	return True	

def printObjAttr(var,obj):
        print var, " = ", vars(obj)[str(var)]
        return var

def printAllObjAttr(obj):
	for item in vars(obj):
		print item, " = ", vars(obj)[str(item)]
	return True	