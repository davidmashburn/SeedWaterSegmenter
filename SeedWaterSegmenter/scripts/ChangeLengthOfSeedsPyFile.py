#!/usr/bin/env python
import os
import wx
import SeedWaterSegmenter as SWS


def forceNewLength(l, newLength, defaultElement=None):
    '''Return a new list of length "newLength" using as many elements from the input list "l" as possible
       If the new length is longer, fill in with elements "defaultElement"'''
    if len(l) >= newLength:
        return l[:newLength]
    elif len(l) < newLength:
        return l + [defaultElement for i in range(newLength - len(l))]


if __name__ == "__main__":
    app = wx.App()

    # Pick a Seeds.py file:
    seedPointsFile = wx.FileSelector("Pick a Seeds.py file")
    if os.path.split(seedPointsFile)[1] != "Seeds.py":
        del seedPointsFile
        wx.MessageBox("You picked a file that was not a Seeds.py !")
        exit()

    # Load all the variables:
    Seeds = SWS.LoadSeedPointsFile(seedPointsFile)
    if Seeds == None:  # Exit if seedList or seedVals is missing
        exit()

    oldLength = len(Seeds.seedList)
    # Choose a new length:
    ndialog = wx.NumberEntryDialog(
        None,
        "Choose a new length for\n"
        + seedPointsFile
        + "\n(old length is "
        + str(oldLength)
        + ")",
        "",
        "Set Seeds.py Length",
        oldLength,
        0,
        2 ** 31 - 1,
    )
    if ndialog.ShowModal() == wx.ID_CANCEL:
        exit()

    newLength = ndialog.GetValue()

    if newLength < oldLength:
        mdialog = wx.MessageDialog(
            None,
            "You are potentially deleting data from\n"
            + seedPointsFile
            + "\n(You might want to make a backup first)"
            + "\nAre you sure you want to do this?",
        )
        if mdialog.ShowModal() == wx.ID_CANCEL:
            exit()

    # Change all the variables:
    seedList = forceNewLength(Seeds.seedList, newLength)
    seedVals = forceNewLength(Seeds.seedVals, newLength)
    walgorithm = (
        forceNewLength(Seeds.walgorithm, newLength, defaultElement="PyMorph")
        if hasattr(Seeds, "walgorithm")
        else ["PyMorph"] * newLength
    )
    woundCenters = (
        forceNewLength(Seeds.woundCenters, newLength)
        if hasattr(Seeds, "walgorithm")
        else [None] * newLength
    )

    # Write back all the variables:
    SWS.WriteSeedPointsFile(
        seedPointsFile, seedList, seedVals, walgorithm, woundCenters
    )
