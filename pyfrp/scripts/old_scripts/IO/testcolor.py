import colorama
colorama.init(autoreset=True)


def printWarning(txt):
	print(colorama.Fore.YELLOW + "WARNING:"), txt

def printError(txt):
	print(colorama.Fore.RED + "ERROR:"), txt

def printGreen(txt):
	print(colorama.Fore.GREEN + txt)
	
printWarning("lalalalal")
printError("bla")
printGreen("Green")
print "lalalalal"