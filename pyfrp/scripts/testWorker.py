#Import modules
from pyfrp.modules import pyfrp_misc_module
from PyQt4 import QtGui, QtCore
import sys

#Simple test function
def testFunc():
	
	#List that is going to be passed along
	l=[3,3,4,56,1,3]
	
	#Call pyfrp function
	print pyfrp_misc_module.remRepeatsList(l)
	
#Main Window
class Window(QtGui.QWidget):

	def __init__(self, parent = None):
	
		QtGui.QWidget.__init__(self, parent)
		
		#Start button
		self.btnStart=QtGui.QPushButton('Start')
		self.btnStart.connect(self.btnStart, QtCore.SIGNAL('clicked()'), self.startPressed)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.btnStart)
		
		self.setLayout(self.hbox)    
		
		self.show()
	
	def startPressed(self):
		
		#Make and start thread
		self.thread = GenericThread()
		self.thread.start()
		
		#Make worker
		self.worker = GenericWorker(testFunc)
		
		#Pass worker to thread
		self.worker.moveToThread(self.thread)
		
		
		self.worker.finished.connect(self.thread.quit)
		
		#Run
		self.worker.start.emit()
		
		
		#self.worker = worker
		
		#thread.wait()
		
		
		

class GenericThread(QtCore.QThread):
	
	def __init__(self):
		QtCore.QThread.__init__(self)
		
	#def __del__(self):
		#self.wait()
		
		
		
				
#Simple worker class        
class GenericWorker(QtCore.QObject):
	
	finished = QtCore.pyqtSignal()
	
	def __init__(self, function, *args, **kwargs):
		super(GenericWorker, self).__init__()

		self.function = function
		self.args = args
		self.kwargs = kwargs
		self.start.connect(self.run)

	start = QtCore.pyqtSignal()

	@QtCore.pyqtSlot()
	def run(self):
		
		self.function(*self.args, **self.kwargs)
		self.finished.emit()
			
#Main
def main():
		    
	app = QtGui.QApplication(sys.argv)
	font=app.font()
	font.setPointSize(12)
	app.setFont(font)
	
	mainWin = Window()
	mainWin.show()
	
	sys.exit(app.exec_())
		
if __name__ == '__main__':
	main()