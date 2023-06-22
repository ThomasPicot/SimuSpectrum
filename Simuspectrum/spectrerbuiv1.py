import matplotlib as mpl
import numpy as np
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from ClassSpectrum.classspectrum import DisplaySpectrum

mpl.use('Qt5Agg')
displayspectrum = DisplaySpectrum()


def openFile(doc):
    with open(doc) as file:
        lis = file.readlines()
        lis = [float(i) * 1e9 for i in lis]
        array = np.array(lis)
        file.close()
        return array


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=0, height=0, dpi=0):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111)
        self.ax.axes.set_ylim(0, 1)
        self.ax.axes.set_xlim(-4.5, 6)
        super(MplCanvas, self).__init__(fig)


class Ui_mainWindow(QtWidgets.QMainWindow):

    def setupUi(self, mainWindow):
        mainWindow.setObjectName("mainWindow")
        mainWindow.resize(655, 849)

        self.centralwidget = QtWidgets.QWidget(mainWindow)
        self.centralwidget.setObjectName("centralwidget")

        self.main_canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.main_canvas.ax.set(xlabel='detuning $\\omega-\\omega_i$ [GHz]', ylabel='Transmission T [su]',
                                title='transmission of rubidium gaz with different cells and temperature ')
        self.toolbar = NavigationToolbar(self.main_canvas, self)

        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(9, 240, 631, 571))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_2.addWidget(self.toolbar)
        self.verticalLayout_2.addWidget(self.main_canvas)

        self.horizontalSlider = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider.setGeometry(QtCore.QRect(240, 50, 201, 21))
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.horizontalSlider.setMinimum(250)
        self.horizontalSlider.setMaximum(500)
        self.horizontalSlider.setSingleStep(1)
        self.horizontalSlider.setValue(300)
        self.horizontalSlider.valueChanged.connect(self.valueUpdateT)

        self.horizontalSlider_2 = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider_2.setGeometry(QtCore.QRect(240, 90, 201, 21))
        self.horizontalSlider_2.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_2.setObjectName("horizontalSlider_2")
        self.horizontalSlider_2.setMinimum(1)
        self.horizontalSlider_2.setMaximum(200)
        self.horizontalSlider_2.setSingleStep(1)
        self.horizontalSlider_2.setValue(75)
        self.horizontalSlider_2.valueChanged.connect(self.valueUpdateLength)

        self.horizontalSlider_3 = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider_3.setGeometry(QtCore.QRect(240, 130, 201, 21))
        self.horizontalSlider_3.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_3.setObjectName("horizontalSlider_3")
        self.horizontalSlider_3.setMinimum(0)
        self.horizontalSlider_3.setMaximum(100)
        self.horizontalSlider_3.setSingleStep(1)
        self.horizontalSlider_3.setValue(28)
        self.horizontalSlider_3.valueChanged.connect(self.valueUpdateFraction)

        self.horizontalSlider_4 = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider_4.setGeometry(QtCore.QRect(20, 210, 201, 21))
        self.horizontalSlider_4.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_4.setObjectName("horizontalSlider_4")
        self.horizontalSlider_4.setMinimum(0)
        self.horizontalSlider_4.setMaximum(100)
        self.horizontalSlider_4.setSingleStep(1)
        self.horizontalSlider_4.setValue(50)
        self.horizontalSlider_4.valueChanged.connect(self.translateExperimentalPlot)

        self.clearButton = QtWidgets.QPushButton(self.centralwidget)
        self.clearButton.setGeometry(QtCore.QRect(540, 120, 81, 21))
        self.clearButton.setObjectName("clearButton")
        self.clearButton.clicked.connect(self.pushClearButton)

        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(90, 50, 131, 21))
        self.label.setObjectName("label")

        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(90, 90, 121, 21))
        self.label_2.setObjectName("label_2")

        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(450, 50, 41, 31))
        self.textBrowser.setObjectName("textBrowser")

        self.textBrowser_2 = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser_2.setGeometry(QtCore.QRect(450, 90, 41, 31))
        self.textBrowser_2.setObjectName("textBrowser_2")

        self.textBrowser_3 = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser_3.setGeometry(QtCore.QRect(450, 130, 41, 31))
        self.textBrowser_3.setObjectName("textBrowser_3")

        self.textBrowser_4 = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser_4.setGeometry(QtCore.QRect(30, 180, 81, 31))
        self.textBrowser_4.setObjectName("textBrowser_4")

        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(90, 130, 131, 21))
        self.label_3.setObjectName("label_3")

        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(250, 210, 171, 21))
        self.label_4.setObjectName("label_4")

        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(250, 210, 151, 21))
        self.label_5.setObjectName("label_5")

        mainWindow.setCentralWidget(self.centralwidget)

        self.menubar = QtWidgets.QMenuBar(mainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuOpen_File = QtWidgets.QMenu(self.menubar)
        self.menuOpen_File.setObjectName("menuOpen_File")
        mainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(mainWindow)
        self.statusbar.setObjectName("statusbar")
        mainWindow.setStatusBar(self.statusbar)

        self.actionOpen_File = QtWidgets.QAction(mainWindow)
        self.actionOpen_File.setObjectName("actionOpen_File")

        self.actionopen_mono_file = QtWidgets.QAction(mainWindow)
        self.actionopen_mono_file.setObjectName("actionopen_mono_file")

        self.menuOpen_File.addAction(self.actionOpen_File)
        self.menubar.addAction(self.menuOpen_File.menuAction())
        self.menuOpen_File.addAction(self.actionopen_mono_file)

        self.actionOpen_File.triggered.connect(self.openFiles)
        self.actionopen_mono_file.triggered.connect(self.openMonoFile)

        self._plot_ref = None
        self.update_ref = None
        self.updateCanvas()

        self.timer = QtCore.QTimer()
        self.timer.setInterval(100)  # in ms
        self.timer.timeout.connect(self.updateCanvas)
        self.timer.start()

        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)

    # returns the slider's value
    def valueUpdateT(self):
        value = str(self.horizontalSlider.value())
        self.textBrowser.setText(value)
        value = float(value)
        return value

    def valueUpdateLength(self):
        value = self.horizontalSlider_2.value()
        value_text = str(value)
        value = float(value / 1000)
        self.textBrowser_2.setText(value_text)
        return value

    def valueUpdateFraction(self):
        value = self.horizontalSlider_3.value()
        value_text = str(value)
        value = float(value)/10
        self.textBrowser_3.setText(value_text)
        return value

    def translateExperimentalPlot(self):
        value = self.horizontalSlider_4.value()
        value = float(value)
        self.textBrowser_4.setText(f"{np.round(-10**-1*((value/10)-5), 2)}GHz")
        return value

    def updateCanvas(self):
        if self._plot_ref is None:
            plot_refs = self.main_canvas.ax.plot(displayspectrum.detuning() * 1e-9,
                                                 displayspectrum.transmission(
                                                     self.valueUpdateFraction(), self.valueUpdateT(),
                                                     self.valueUpdateLength()), label='theory')
            self.main_canvas.ax.legend(loc="center right")
            self._plot_ref = plot_refs[0]
        else:
            self._plot_ref.set_ydata(displayspectrum.transmission(self.valueUpdateFraction(),
                                                                  self.valueUpdateT(), self.valueUpdateLength()))

        self.main_canvas.draw()

    # you just need to return the file's name
    def openFiles(self):
        txtin = str(QFileDialog.getOpenFileNames())
        if txtin == str(([], '')):
            pass
        else:
            txtout = str(QFileDialog.getOpenFileNames())
            if txtout == str(([], '')):
                pass
            else:
                n = len(txtin)
                txtin = txtin[3:n - 20]
                n = len(txtout)
                txtout = txtout[3:n - 20]
                inpute = openFile(txtin)
                output = openFile(txtout)
                transm = np.array([i / o for o, i in zip(output, inpute)])
                x = openFile("detuning.txt")
                self.main_canvas.ax.plot(x / 1e9, transm / transm[0], label='experiment')
                self.main_canvas.ax.legend(loc="center right")

    def openMonoFile(self):
        experimental_spectrum = str(QFileDialog.getOpenFileNames())
        if experimental_spectrum == str(([], '')):
            pass
        else:
            n = len(experimental_spectrum)
            experimental_spectrum = openFile(experimental_spectrum[3:n - 20])
            x = openFile("detuning.txt")
            self.main_canvas.ax.plot((x/1e9)-(self.translateExperimentalPlot()/100)+0.5,
                                     experimental_spectrum / experimental_spectrum[0], label='experimental')
            self.main_canvas.ax.legend()

    def pushClearButton(self):
        self.setupUi(mainWindow=mainWindow)

    def retranslateUi(self, mainWindow):
        _translate = QtCore.QCoreApplication.translate
        mainWindow.setWindowTitle(_translate("mainWindow", "Rb87 Transmission "))
        self.clearButton.setText(_translate("MainWindow", "clear button"))
        self.label.setText(_translate("mainWindow", "Temperature [K]"))
        self.label_4.setText(_translate("mainWindow", "translate experimental plot"))
        self.label_2.setText(_translate("mainWindow", "Length L [mm]"))
        self.label_3.setText(_translate("mainWindow", "fraction Rb 87/85 [%]"))
        self.label_4.setText(_translate("mainWindow", "translation in GHz"))
        self.menuOpen_File.setTitle(_translate("mainWindow", "File"))
        self.actionOpen_File.setText(_translate("mainWindow", "select out/input Files "))
        self.actionopen_mono_file.setText(_translate("mainWindow", "open mono-file"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()
    ui = Ui_mainWindow()
    ui.setupUi(mainWindow)
    mainWindow.show()
    sys.exit(app.exec_())


