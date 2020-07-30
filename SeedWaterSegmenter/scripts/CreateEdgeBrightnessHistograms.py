#!/usr/bin/env python
"""Plot histograms of the edgeBrightnessValues
   which is generated with the RunCalculations button.
   By default, this only plots histograms for
   the first and last frames."""
import matplotlib.pyplot as plt
import wx

app = wx.App()

ebFile = wx.FileSelector(
    "Select an edgeBrightnessValues.csv...", default_extension=".csv", wildcard="*.csv"
)
brightnessValues = [map(int, i.split(",")) for i in open(ebFile)]

plt.ioff()
plt.clf()
plt.hist(brightnessValues[0], bins=20)
plt.hist(brightnessValues[-1], bins=20)
plt.show()
