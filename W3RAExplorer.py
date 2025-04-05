from qgis.utils import iface
from qgis.gui import QgsMapTool
from qgis.core import QgsProject
from PyQt5.QtWidgets import QAction, QFileDialog, QInputDialog, QMessageBox, QComboBox, QDialog, QVBoxLayout, QPushButton, QLineEdit, QLabel, QHBoxLayout
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas  # updated: for embedded plots
from matplotlib.figure import Figure
import numpy as np
import os
import re
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d  # added: for Gaussian smoothing
from qgis.gui import QgsMapTool, QgsVertexMarker  # <-- added QgsVertexMarker
from PyQt5.QtCore import Qt  # <-- needed for Qt.red


class InteractivePlotDialog(QDialog):  # updated: persistent regression interface
    def __init__(self, x, y, anomalies, var_name, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Time Series Analysis")
        self.resize(800, 600)

        self.x = x
        self.y = y
        self.anomalies = anomalies
        self.var_name = var_name

        self.layout = QVBoxLayout(self)

        self.combo = QComboBox()
        self.combo.addItems(["None", "Linear", "Polynomial", "Exponential", "Gaussian Smoothing", "Fourier"])
        self.combo.currentIndexChanged.connect(self.update_plot)

        self.param_input = QLineEdit()
        self.param_input.setPlaceholderText("Degree (Poly), Sigma (Gaussian), Terms (Fourier)")
        self.param_input.textChanged.connect(self.update_plot)

        self.layout.addWidget(QLabel("Regression Type:"))
        self.layout.addWidget(self.combo)
        self.layout.addWidget(self.param_input)

        self.fig = Figure(figsize=(10, 6))
        self.canvas = FigureCanvas(self.fig)
        self.layout.addWidget(self.canvas)

        self.plot_data()

    def plot_data(self):
        self.fig.clear()
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)

        ax1.plot(self.x, self.y, 'o-', label=self.var_name, color='red', markersize=4)
        ax1.set_title(f"{self.var_name} Time Series")
        ax1.set_ylabel("Value")
        ax1.grid(True)

        ax2.plot(self.x, self.anomalies, 'o--', label="Anomalies", color='blue', markersize=4)
        ax2.set_title("Anomalies")
        ax2.set_xlabel("Time (days)")
        ax2.set_ylabel("Anomaly")
        ax2.grid(True)

        self.ax1 = ax1
        self.canvas.draw()
        self.update_plot()

    def update_plot(self):
        reg_type = self.combo.currentText()
        param = self.param_input.text()

        while len(self.ax1.lines) > 1:
          line = self.ax1.lines[-1]
          line.remove()
      # self.ax1.lines = self.ax1.lines[:1]  # keep only raw signal
        x = self.x
        y = self.y

        try:
            if reg_type == "Linear":
                coeffs = np.polyfit(x, y, 1)
                y_fit = np.polyval(coeffs, x)
                self.ax1.plot(x, y_fit, label="Linear Fit", linestyle='--')

            elif reg_type == "Polynomial":
                deg = int(param) if param else 2
                coeffs = np.polyfit(x, y, deg)
                y_fit = np.polyval(coeffs, x)
                self.ax1.plot(x, y_fit, label=f"Poly deg={deg}", linestyle='--')

            elif reg_type == "Exponential":
                def model(x, a, b, c): return a * np.exp(b * x) + c
                popt, _ = curve_fit(model, x, y, maxfev=10000)
                self.ax1.plot(x, model(x, *popt), label="Exp Fit", linestyle='--')

            elif reg_type == "Gaussian Smoothing":
                sigma = float(param) if param else 3
                y_smooth = gaussian_filter1d(y, sigma)
                self.ax1.plot(x, y_smooth, label=f"Gaussian σ={sigma}", linestyle='--')

            elif reg_type == "Fourier":
                terms = int(param) if param else 5
                fft_coeffs = np.fft.rfft(y)
                fft_coeffs[terms:] = 0
                y_ifft = np.fft.irfft(fft_coeffs)
                self.ax1.plot(x, y_ifft, label=f"Fourier ({terms} terms)", linestyle='--')

        except Exception as e:
            print("Regression update error:", e)

        self.ax1.legend()
        self.canvas.draw()

class PointClickTool(QgsMapTool):
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        self.marker = None

    def canvasReleaseEvent(self, event):
        point = self.canvas.getCoordinateTransform().toMapCoordinates(event.pos())
        lon, lat = point.x(), point.y()

        # Remove any existing marker
        if self.marker:
            self.canvas.scene().removeItem(self.marker)

        # Add a new marker at the clicked location
        self.marker = QgsVertexMarker(self.canvas)
        self.marker.setCenter(point)
        self.marker.setColor(Qt.red)
        self.marker.setIconSize(12)
        self.marker.setIconType(QgsVertexMarker.ICON_CROSS)
        self.marker.setPenWidth(3)

        input_file = QgsProject.instance().readEntry("W3RAExplorer", "NetCDF_Path")[0]
        if not input_file or not os.path.exists(input_file):
            QMessageBox.critical(None, "Error", "NetCDF file not found. Please load a valid file.")
            return

        try:
            dataset = nc.Dataset(input_file, "r")
            pattern = re.compile(r"(.+)_\d{4}$")
            base_vars = sorted(set(
                pattern.match(v).group(1)
                for v in dataset.variables
                if pattern.match(v) and dataset.variables[v].ndim == 3
            ))

            if not base_vars:
                QMessageBox.warning(None, "No Variables", "No time-dependent 3D variables found.")
                return

            var_name, ok = QInputDialog.getItem(None, "Select Variable", "Variable:", base_vars, 0, False)
            if not ok:
                return

            years = sorted(
                int(v.split("_")[-1])
                for v in dataset.variables
                if v.startswith(var_name + "_") and dataset.variables[v].ndim == 3
            )

            lat_arr = dataset.variables["lat"][:]
            lon_arr = dataset.variables["lon"][:]
            lat_idx = (np.abs(lat_arr - lat)).argmin()
            lon_idx = (np.abs(lon_arr - lon)).argmin()

            all_times = []
            all_values = []

            for year in years:
                vname = f"{var_name}_{year}"
                tname = f"time_{year}"
                if vname in dataset.variables and tname in dataset.variables:
                    times = dataset.variables[tname][:]
                    values = dataset.variables[vname][:, lat_idx, lon_idx]
                    all_times.extend(times + (year - years[0]) * 365)
                    all_values.extend(values)

            dataset.close()

            if not all_values:
                QMessageBox.warning(None, "No Data", "No data found at this location.")
                return

            x = np.array(all_times)
            y = np.array(all_values)
            anomalies = y - np.mean(y)

            dlg = InteractivePlotDialog(x, y, anomalies, var_name)
            dlg.exec_()

        except Exception as e:
            QMessageBox.critical(None, "Error", str(e))

class W3RAExplorer:
    def __init__(self, iface):
        self.iface = iface
        self.tool = None
        self.action = None

    def initGui(self):
        self.tool = PointClickTool(self.iface.mapCanvas())
        self.action = QAction("W3RA Explorer", self.iface.mainWindow())
        self.action.triggered.connect(self.activate_plugin)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&W3RA Explorer", self.action)

    def activate_plugin(self):
        input_file, _ = QFileDialog.getOpenFileName(None, "Select NetCDF File", "", "NetCDF (*.nc)")
        if input_file:
            QgsProject.instance().writeEntry("W3RAExplorer", "NetCDF_Path", input_file)
            self.iface.mapCanvas().setMapTool(self.tool)

    def unload(self):
        if self.action:
            self.iface.removeToolBarIcon(self.action)
            self.iface.removePluginMenu("&W3RA Explorer", self.action)
        if self.tool:
            self.iface.mapCanvas().unsetMapTool(self.tool)
