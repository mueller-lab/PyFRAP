#Script to test how class inheritance actually works

"""
    ---first---
    |          |
    |          |
  second     third
    |          |
    |          |
    ---fourth---
"""

#Define classes
class first(object):

	def __init__(self,name):
		self.name=name
		print name
		
	def sayFirst(self):
		print "first"
		
class second(first):

	def __init__(self,name,typ):
		
		first.__init__(self,name)
		self.typ=typ
		
	def saySecond(self):
		print "second"	
		
class third(first):

	def __init__(self,name,size):
		first.__init__(self,name)
		self.size=size
	
	def sayThird(self):
		print "third"	
		
class fourth(second,third):

	def __init__(self,name,cool,size,typ):
		
		second.__init__(self,name,typ)
		third.__init__(self,name,size)
		
		self.cool=cool
	
	def sayFourth(self):
		print "fourth"

	#Overwrite function
	def sayFirst(self):
		print "No"

#Run script

print 
print 



x=fourth("fromFourth","definitelycool",2.5,1)

print "====="
print x.name
print x.size
print x.typ
print x.cool
print

x.sayFirst()
x.saySecond()
x.sayThird()
x.sayFourth()




				

		