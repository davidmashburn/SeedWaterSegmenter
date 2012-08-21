#! python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
#from wx_py import __file__

import imp

name = 'SeedWaterSegmenter'

modulePath = imp.find_module(name)[1]

if sys.platform == "win32":
    from postinstall import get_special_folder_path
    DESKTOP_FOLDER = get_special_folder_path("CSIDL_DESKTOPDIRECTORY")

    if sys.argv[1] == '-install':
        lnkName = name+'.lnk'
        create_shortcut(
            os.path.join(sys.prefix, 'pythonw.exe'), # program
            name, # description
            lnkName, # filename
            os.path.join(modulePath, shellName+'.py'), # parameters
            '', # workdir
            os.path.join(modulePath,'icons',name+'.ico'), # iconpath
        )
        # move shortcut from current directory to DESKTOP_FOLDER
        shutil.move(os.path.join(os.getcwd(), lnkName),
                    os.path.join(DESKTOP_FOLDER, lnkName))
        # tell windows installer that we created another
        # file which should be deleted on uninstallation
        file_created(os.path.join(DESKTOP_FOLDER, lnkName))

    if sys.argv[1] == '-remove':
        pass
        # This will be run on uninstallation. Nothing to do.

