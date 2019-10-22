import sys
import numpy as np
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.Qt import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MOT(QMainWindow):
	def __init__(self):
		super().__init__()
		self.width = 640
		self.height = 400
		self.top = 50
		self.left = 50
		self.title = 'MOT'
		self.initUI()

	def initUI(self):
		self.setWindowTitle(self.title)
		self.setGeometry(self.left,self.top,self.width,self.height)

		self.plot_land = PlotWorld(self,width=5,height=4)
		self.plot_land.move(0,0)

		self.tools = QDockWidget(self)
		self.toolgroup = QGroupBox(self)
		layout = QGridLayout()

		self.num_lab = QLabel('# of Coils\nper Set')
		self.num_val = QLineEdit(str(self.plot_land.N))
		layout.addWidget(self.num_lab,0,1)
		layout.addWidget(self.num_val,0,2)
		self.cur_lab = QLabel('Current(A)')
		self.cur_val = QLineEdit(str(self.plot_land.I))
		layout.addWidget(self.cur_lab,1,1)
		layout.addWidget(self.cur_val,1,2)
		self.sep_lab = QLabel('Major\nOffset(cm)')
		self.sep_val = QLineEdit(str(self.plot_land.A))
		layout.addWidget(self.sep_lab,2,1)
		layout.addWidget(self.sep_val,2,2)
		self.sep2_lab = QLabel('Minor\nOffset(cm)')
		self.sep2_val = QLineEdit(str(self.plot_land.Y))
		layout.addWidget(self.sep2_lab,3,1)
		layout.addWidget(self.sep2_val,3,2)
		self.rmin_lab = QLabel('rmin(cm)')
		self.rmin_val = QLineEdit(str(self.plot_land.rmin))
		layout.addWidget(self.rmin_lab,4,1)
		layout.addWidget(self.rmin_val,4,2)
		self.rmax_lab = QLabel('rmax(cm)')
		self.rmax_val = QLineEdit(str(self.plot_land.rmax))
		layout.addWidget(self.rmax_lab,5,1)
		layout.addWidget(self.rmax_val,5,2)
		self.plot_button = QPushButton('PLOT')
		layout.addWidget(self.plot_button,6,1,1,2)
		self.plot_button.clicked.connect(self.plot_butt_clicked)

		self.grad_lab = QLabel('Gradient at\nz = 0 (Gauss/cm)')
		self.grad_val = QLabel(str(self.plot_land.Gzero))
		layout.addWidget(self.grad_lab,7,1)
		layout.addWidget(self.grad_val,7,2)

		self.blank_lab = QLabel(' '*180)
		layout.addWidget(self.blank_lab,0,0)

		self.toolgroup.setLayout(layout)
		self.tools.setWidget(self.toolgroup)
		self.addDockWidget(Qt.RightDockWidgetArea,self.tools)

		self.num_val.returnPressed.connect(self.plot_button.click)
		self.cur_val.returnPressed.connect(self.plot_button.click)
		self.sep_val.returnPressed.connect(self.plot_button.click)
		self.sep2_val.returnPressed.connect(self.plot_button.click)
		self.rmin_val.returnPressed.connect(self.plot_button.click)
		self.rmax_val.returnPressed.connect(self.plot_button.click)

		self.show()

	def plot_butt_clicked(self):
		self.plot_land.ax1.clear()
		self.plot_land.ax2.clear()
		self.plot_land.N = int(self.num_val.text())
		self.plot_land.I = float(self.cur_val.text())
		self.plot_land.A = float(self.sep_val.text())
		self.plot_land.Y = float(self.sep2_val.text())
		self.plot_land.rmin = float(self.rmin_val.text())
		self.plot_land.rmax = float(self.rmax_val.text())
		self.plot_land.plot()
		self.grad_val.setText(str(self.plot_land.Gzero))
		self.update()



class PlotWorld(FigureCanvas):
	def __init__(self,parent=None,width=5,height=4,dpi=100):
		fig = Figure(figsize=(width,height),dpi=dpi)
		self.axes = fig.add_subplot(211)
		FigureCanvas.__init__(self,fig)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)
		self.n = 1000 # resolution of z axis
		self.uo = 1e-7 * 100 * 1e4 # (Gauss*cm/A) vacuum permeability / 4*pi

		self.N = 6 # Number of Coils
		self.I = 6 # (A) Current Through Coils
		self.Z = np.linspace(-1,1,self.n) # (cm) position on z axis (trap at z = 0)
		self.A = .655 # (cm) distance from coil center to z = 0
		self.Y = .091 # (cm) 1/2 distance between two sides of the coil

		self.rmax = 1.4605 # (cm) maximum radius to fit into chamber (lower than 3.81)
		self.rmin = .6985 # (cm) minimum radius that can easily be machined and glued 
		
		self.plot()

	def plot(self):
		self.R = np.linspace(self.rmin,self.rmax,num=self.N) # (cm) radius for each coil
		self.B = np.zeros((self.N,self.n))
		self.Btot = np.zeros(self.n)
		self.Gtot = np.zeros(self.n)
		self.Gzero = 0

		for r in range(len(self.R)):
			for z in range(len(self.Z)):
				self.B[r,z] = (
					self.uo*2*np.pi*self.I*self.R[r]**2/(((self.Z[z]-self.A-self.Y)**2+self.R[r]**2)**(3/2)) -
				  	self.uo*2*np.pi*self.I*self.R[r]**2/(((self.Z[z]+self.A-self.Y)**2+self.R[r]**2)**(3/2)) +
				  	self.uo*2*np.pi*self.I*self.R[r]**2/(((self.Z[z]-self.A+self.Y)**2+self.R[r]**2)**(3/2)) - 
				  	self.uo*2*np.pi*self.I*self.R[r]**2/(((self.Z[z]+self.A+self.Y)**2+self.R[r]**2)**(3/2))
				  	)
		self.Btot = self.B.sum(axis=0) # (Gauss) Total Magnetic Field along z axis
		self.Gtot = np.gradient(self.Btot)*self.n/2 # (Gauss/cm) Total Magnetic Field GRADIENT along z axis
		self.Gzero = self.Gtot[len(self.Gtot)//2]

		ymin1 = min(self.Btot)
		ymax1 = max(self.Btot)
		ymin2 = min(self.Gtot)
		ymax2 = max(self.Gtot)
		coil_lines = [self.A+self.Y,self.A-self.Y,-self.A+self.Y,-self.A-self.Y]
		colors = ['red','red','blue','blue']
		self.ax1 = self.figure.add_subplot(211)
		self.ax2 = self.figure.add_subplot(212)
		self.ax1.plot(self.Z,self.Btot)
		self.ax1.vlines(coil_lines,ymin1,ymax1,colors=colors)
		self.ax2.plot(self.Z,self.Gtot)
		self.ax2.vlines(coil_lines,ymin2,ymax2,colors=colors)
		self.draw()


if __name__ == '__main__':
	app = QApplication(sys.argv)
	ex = MOT()
	sys.exit(app.exec_())