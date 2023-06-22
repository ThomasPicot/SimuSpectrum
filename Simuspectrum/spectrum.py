# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui/spectrum.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from pyqtgraph import PlotWidget
import numpy as np
from scipy.signal import find_peaks

from Classes.classspectrum import SpectrumRubidiumD2Line
from Classes.ClassLaser import USBScope
ds = SpectrumRubidiumD2Line()


class Ui_SpectrumView(object):

    def __init__(self, addr=None):
        self.addr = addr
        self.timer = QtCore.QTimer()
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.plotTheory)
        self.timer.start()

        if self.addr is not None:
            self.scope = USBScope(addr=self.addr)
        else:
            print("no scope selected")

    def setupUi(self, SpectrumView):
        SpectrumView.setObjectName("SpectrumView")
        SpectrumView.resize(1247, 986)
        self.centralwidget = QtWidgets.QWidget(SpectrumView)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.main_layout = QtWidgets.QGridLayout()
        self.main_layout.setContentsMargins(1, 1, 1, 1)
        self.main_layout.setObjectName("main_layout")
        self.x_translation = QtWidgets.QSlider(self.centralwidget)
        self.x_translation.setOrientation(QtCore.Qt.Horizontal)
        self.x_translation.setObjectName("x_translation")
        self.main_layout.addWidget(self.x_translation, 2, 0, 1, 1)
        self.sliders_layout = QtWidgets.QGridLayout()
        self.sliders_layout.setObjectName("sliders_layout")

        self.mixture_slider = QtWidgets.QSlider(self.centralwidget)
        self.mixture_slider.setOrientation(QtCore.Qt.Horizontal)
        self.mixture_slider.setObjectName("mixture_slider")
        self.mixture_slider.setMinimum(0)
        self.mixture_slider.setMaximum(1000)
        self.mixture_slider.setSingleStep(1)
        self.mixture_slider.setValue(280)
        self.mixture_slider.valueChanged.connect(self.valueUpdateFraction)

        self.sliders_layout.addWidget(self.mixture_slider, 2, 1, 1, 1)
        self.lineEdit_temp = QtWidgets.QLineEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_temp.sizePolicy().hasHeightForWidth())
        self.lineEdit_temp.setSizePolicy(sizePolicy)
        self.lineEdit_temp.setObjectName("lineEdit_temp")
        self.sliders_layout.addWidget(self.lineEdit_temp, 0, 0, 1, 1)
        self.lineEdit_mixture = QtWidgets.QLineEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_mixture.sizePolicy().hasHeightForWidth())
        self.lineEdit_mixture.setSizePolicy(sizePolicy)
        self.lineEdit_mixture.setObjectName("lineEdit_mixture")
        self.sliders_layout.addWidget(self.lineEdit_mixture, 2, 0, 1, 1)
        self.lineEdit_length = QtWidgets.QLineEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_length.sizePolicy().hasHeightForWidth())
        self.lineEdit_length.setSizePolicy(sizePolicy)
        self.lineEdit_length.setObjectName("lineEdit_length")
        self.sliders_layout.addWidget(self.lineEdit_length, 1, 0, 1, 1)

        self.length_slider = QtWidgets.QSlider(self.centralwidget)
        self.length_slider.setOrientation(QtCore.Qt.Horizontal)
        self.length_slider.setObjectName("length_slider")
        self.length_slider.setMinimum(0)
        self.length_slider.setMaximum(200)
        self.length_slider.setSingleStep(1)
        self.length_slider.setValue(75)
        self.length_slider.valueChanged.connect(self.valueUpdateLength)

        self.sliders_layout.addWidget(self.length_slider, 1, 1, 1, 1)

        self.temp_slider = QtWidgets.QSlider(self.centralwidget)
        self.temp_slider.setOrientation(QtCore.Qt.Horizontal)
        self.temp_slider.setObjectName("temp_slider")
        self.temp_slider.setMaximum(500)
        self.temp_slider.setMinimum(300)
        self.temp_slider.setSingleStep(1)
        self.temp_slider.setValue(300)
        self.temp_slider.valueChanged.connect(self.valueUpdateTemp)

        self.sliders_layout.addWidget(self.temp_slider, 0, 1, 1, 1)
        self.temp_label = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.temp_label.sizePolicy().hasHeightForWidth())
        self.temp_label.setSizePolicy(sizePolicy)
        self.temp_label.setObjectName("temp_label")
        self.sliders_layout.addWidget(self.temp_label, 0, 2, 1, 1)
        self.length_label = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.length_label.sizePolicy().hasHeightForWidth())
        self.length_label.setSizePolicy(sizePolicy)
        self.length_label.setObjectName("length_label")
        self.sliders_layout.addWidget(self.length_label, 1, 2, 1, 1)
        self.mixture_label = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mixture_label.sizePolicy().hasHeightForWidth())
        self.mixture_label.setSizePolicy(sizePolicy)
        self.mixture_label.setObjectName("mixture_label")
        self.sliders_layout.addWidget(self.mixture_label, 2, 2, 1, 1)
        self.main_layout.addLayout(self.sliders_layout, 0, 0, 1, 1)
        self.display_widget = PlotWidget(self.centralwidget)
        self.display_widget.setObjectName("display_widget")
        self.display_widget.setXRange(-4, 6)
        self.display_widget.setYRange(0, 1)
        self.display_widget.setLabel(axis='left', text='Transmission')
        self.display_widget.setLabel(axis='bottom', text='detuning [GHz]')
        self.main_layout.addWidget(self.display_widget, 3, 0, 1, 1)
        self.label_translation = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_translation.sizePolicy().hasHeightForWidth())
        self.label_translation.setSizePolicy(sizePolicy)
        self.label_translation.setObjectName("label_translation")
        self.main_layout.addWidget(self.label_translation, 1, 0, 1, 1)
        self.gridLayout_2.addLayout(self.main_layout, 1, 0, 1, 1)
        self.clear_button = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clear_button.sizePolicy().hasHeightForWidth())
        self.clear_button.setSizePolicy(sizePolicy)
        self.clear_button.setObjectName("pushButton")
        self.gridLayout_2.addWidget(self.clear_button, 0, 0, 1, 1)
        SpectrumView.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(SpectrumView)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1247, 22))
        self.menubar.setObjectName("menubar")
        self.menuopen_File = QtWidgets.QMenu(self.menubar)
        self.menuopen_File.setObjectName("menuopen_File")
        SpectrumView.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(SpectrumView)
        self.statusbar.setObjectName("statusbar")
        SpectrumView.setStatusBar(self.statusbar)
        self.actioninput = QtWidgets.QAction(SpectrumView)
        self.actioninput.setObjectName("actioninput")
        self.actionoutpout = QtWidgets.QAction(SpectrumView)
        self.actionoutpout.setObjectName("actionoutpout")
        self.menubar.addAction(self.menuopen_File.menuAction())

        self.retranslateUi(SpectrumView)
        QtCore.QMetaObject.connectSlotsByName(SpectrumView)
        self.plotTheory()

    def retranslateUi(self, SpectrumView):
        _translate = QtCore.QCoreApplication.translate
        SpectrumView.setWindowTitle(_translate("SpectrumView", "Spectrum"))
        self.temp_label.setText(_translate("SpectrumView", "Temperature [K]"))
        self.length_label.setText(_translate("SpectrumView", "lentgth [mm]"))
        self.mixture_label.setText(_translate("SpectrumView", "Fraction Rb 87/85"))
        self.label_translation.setText(_translate("SpectrumView", "Translation experimental graph"))
        self.clear_button.setText(_translate("SpectrumView", "Clear Button"))
        self.menuopen_File.setTitle(_translate("SpectrumView", "open File"))
        self.actioninput.setText(_translate("SpectrumView", "input"))
        self.actionoutpout.setText(_translate("SpectrumView", "outpout"))

    def valueUpdateTemp(self):
        value = str(self.temp_slider.value())
        self.lineEdit_temp.setText(value)
        value = float(value)
        return value

    def valueUpdateFraction(self):
        value = self.mixture_slider.value()
        value_text = str(float(value)/10)
        value = float(value)/10
        self.lineEdit_mixture.setText(value_text)
        return value

    def valueUpdateLength(self):
        value = self.length_slider.value()
        value_text = str(value)
        value = float(value / 1000)
        self.lineEdit_length.setText(value_text)
        return value

    def plotTheory(self):
        self.display_widget.plot(ds.detuning() * 1e-9, ds.transmission(self.valueUpdateFraction(),
                                                                        self.valueUpdateTemp(), self.valueUpdateLength())
                                 , clear=True)

    def plotExpSpectrum(self, channel=1):
        transmission = self.scope.get_waveform(channels=channel)[0]
        transmission = transmission/np.amax(transmission)
        peaks = find_peaks(transmission, distance=3000, width=1000)

    def get_data_with_thread(self):
        """
        call the method that executes the thread.
        :return:
        """
        self.worker.get_data()


class WorkerThread(QtCore.QThread):
    output = QtCore.pyqtSignal(np.ndarray)

    def __init__(self, fun, parent=None):
        """

        :param fun: enter the function you want to let run in another thread, in our case it is the method that plots
        in time the saturated susceptibility cell.
        :param parent: None
        """
        QtCore.QThread.__init__(self, parent)
        self.fun = fun

    def get_data(self):
        """
        call the run (the start is an attribute from QThread)
        :return:
        """
        try:
            self.start()
        except:
            pass

    def run(self):
        """
        run the function. must never be called
        :return: None
        """
        try:
            self.fun()
            self.output.emit(self.fun)
        except:
            pass


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SpectrumView = QtWidgets.QMainWindow()
    ui = Ui_SpectrumView()
    ui.setupUi(SpectrumView)
    SpectrumView.show()
    sys.exit(app.exec_())
