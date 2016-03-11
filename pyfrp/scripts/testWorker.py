#Import modules
from pyfrp.modules import pyfrp_misc_module
from PyQt4 import QtGui, QtCore
import sys
import time

#Simple test function
def testFunc(sig):
	
	#List that is going to be passed along
	l=[3,3,4,56,1,3]
	
	#Call pyfrp function
	print pyfrp_misc_module.remRepeatsList(l)
	
	for i in range(100):
		time.sleep(0.03)	
		sig.emit(i)
	
#Main Window
class Window(QtGui.QWidget):

	def __init__(self, parent = None):
	
		QtGui.QWidget.__init__(self, parent)
		
		#Start button
		self.btnStart=QtGui.QPushButton('Start')
		self.btnStart.connect(self.btnStart, QtCore.SIGNAL('clicked()'), self.startPressed)
		
		#Progressbar
		self.progressbar = QtGui.QProgressBar()
		self.progressbar.setMinimum(1)
		self.progressbar.setMaximum(100)
		
		#Layout
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.btnStart)
		self.hbox.addWidget(self.progressbar)
		
		self.setLayout(self.hbox)    
		
		self.show()
	
	def startPressed(self):
		
		#Make and start thread
		self.thread = GenericThread()
		self.thread.start()
		
		#Make worker
		self.worker = GenericWorker(testFunc,self.thread.progressSignal)
		
		#Pass worker to thread
		self.worker.moveToThread(self.thread)
		
		#Connect
		self.worker.finished.connect(self.quitThread)
		self.thread.progressSignal.connect(self.updateProgressDialog)
		
		#Run
		self.worker.start.emit()
		
	def quitThread(self):
		self.thread.quit()
	
	def updateProgressDialog(self,n):
		self.progressbar.setValue(n)
	
	#def closeEvent(self, event):
			
		#self.thread.__del__()	
			
		#self.worker.__del__()
		
		#return
		
class GenericThread(QtCore.QThread):
	
	progressSignal = QtCore.pyqtSignal(int)
	
	def __init__(self):
		QtCore.QThread.__init__(self)
		
	def __del__(self):
		self.wait()
		
		
		
				
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