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
    try:
        # When running a binary installer, special functions like
        # get_special_folder_path, create_shortcut, and file_created
        # are "magically" available; see
        # http://docs.python.org/2/distutils/builtdist.html
        DESKTOP_FOLDER = get_special_folder_path("CSIDL_DESKTOPDIRECTORY")
    except:
        # but if we are not running it through the binary installer, we need to make calls to PyWin instead
        import win32com.client
        from win32com.shell import shell, shellcon
        
        DESKTOP_FOLDER = shell.SHGetFolderPath( 0,
                                                shellcon.CSIDL_DESKTOPDIRECTORY,
                                                None,0 )
        
        def create_shortcut( target, description, filename,
                             arguments, workdir, iconpath ):
            '''Make a shortcut with direct calls to Windows'''
            shell = win32com.client.Dispatch('WScript.Shell')
            shortcut = shell.CreateShortCut(filename)
            shortcut.TargetPath = target
            shortcut.Description = description
            shortcut.Arguments = arguments
            shortcut.WorkingDirectory = workdir
            shortcut.IconLocation = iconpath
            #shortcut.FullName
            # shortcut.Hotkey
            # shortcut.WindowStyle
            shortcut.Save()
        def file_created(filename):
            '''This is only used for uninstallation, and this will never run during installation'''
            pass
    
    if len(sys.argv)<2:
        print 'No option specified. Please append either "-install" or "-remove"'
    else:
        if sys.argv[1] == '-install':
            lnkName = name+'.lnk'
            pythonw = os.path.join( os.path.split(sys.executable)[0], 'pythonw.exe')
            create_shortcut(
                os.path.join(pythonw), # target
                name, # description
                lnkName, # filename
                os.path.join(modulePath, name+'.py'), # arguments
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
