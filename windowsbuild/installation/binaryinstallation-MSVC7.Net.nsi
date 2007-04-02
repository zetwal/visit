; Script generated by the HM NIS Edit Script Wizard.

##############################################################################
#
# File: binaryinstallation-MSVC7.Net.nsi
#
# Purpose: This file contains the instructions that the NSIS installer needs
#          in order to create an installation program for VisIt.
#
# Programmer: Brad Whitlock
# Date:       Wed Aug 18 15:46:47 PST 2004
#
# Modifications:
#  Brad Whitlock, Thu Sep 23 09:39:53 PDT 2004
#  Updated for 1.3.5.
#
#  Brad Whitlock, Wed Nov 3 14:07:03 PST 2004
#  Updated for 1.4.
#
#  Brad Whitlock, Wed Jan 5 17:24:18 PST 2005
#  Updated for 1.4.1.
#
#  Brad Whitlock, Thu Feb 24 16:09:30 PST 2005
#  Updated for 1.4.2. I also added more configuration screens that allow
#  the user to pick a default database format.
#
#  Brad Whitlock, Tue May 10 14:13:33 PST 2005
#  Updated for 1.4.3.
#
#  Brad Whitlock, Mon Jun 6 17:01:48 PST 2005
#  Added support for setting VisIt's install location as a Java preference
#  using our VIkit plugin.
#
#  Brad Whitlock, Fri Jun 24 10:35:49 PDT 2005
#  Added support for saving movies at the resolution stored in the
#  session file.
#
#  Brad Whitlock, Fri Jul 8 15:40:16 PST 2005
#  Added code to move config files if there are any from older versions
#  of VisIt.
#
#  Brad Whitlock, Thu Aug 18 12:00:55 PDT 2005
#  Updated for version 1.4.5.
#
#  Brad Whitlock, Tue Nov 22 13:53:11 PST 2005
#  Updated for 1.5.
#
#  Brad Whitlock, Wed Feb 1 09:33:46 PDT 2006
#  Updated for 1.5.1
#
#  Brad Whitlock, Wed Mar 15 9:57:34 PDT 2006
#  Updated for 1.5.2. I also added support for NERSC and ORNL network
#  configs as well as a new screen that lets you pick your default
#  bank for parallel jobs.
#
#  Brad Whitlock, Tue Jun 13 15:12:56 PST 2006
#  Updated for 1.5.3.
#
#  Kathleen Bonnell, Wed Oct 18 09:11:31 PDT 2006
#  Updated for 1.5.4
#
#  Brad Whitlock, Tue Nov 21 16:41:49 PST 2006
#  Upped the version number and added movietemplates directories
#
##############################################################################

; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "VisIt"
!define PRODUCT_VERSION "1.5.5"
!define PRODUCT_PUBLISHER "LLNL"
!define PRODUCT_WEB_SITE "http://www.llnl.gov/visit"
!define PRODUCT_DIR_REGKEY "Software\Microsoft\Windows\CurrentVersion\App Paths\visit${PRODUCT_VERSION}.exe"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME} ${PRODUCT_VERSION}"
!define PRODUCT_UNINST_ROOT_KEY "HKLM"

# Define macros for Qt.
!define QTPATH "C:\Qt\3.0.2"
!define QTDLL  "qt-mt302.dll"

SetCompressor bzip2

; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "..\resources\visit.ico"
!define MUI_UNICON "..\resources\visit.ico"
#!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\modern-uninstall.ico"

ReserveFile "NetworkConfig.ini"
ReserveFile "WantDatabasePlugin.ini"
ReserveFile "DefaultDatabasePlugin.ini"
ReserveFile "ChooseInstallDevelopmentFiles.ini"
ReserveFile "ClickInstall.ini"

; Reserve files
!insertmacro MUI_RESERVEFILE_INSTALLOPTIONS

; Welcome page
!insertmacro MUI_PAGE_WELCOME
; License page
!insertmacro MUI_PAGE_LICENSE "copyright.txt"
; Directory page
!insertmacro MUI_PAGE_DIRECTORY
; Custom
page custom ChooseNetworkConfig
page custom ChooseParallelBank
page custom WantDefaultDatabasePlugin
page custom ChooseDefaultDatabasePlugin
page custom ChooseInstallDevelopmentFiles

; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
!define MUI_FINISHPAGE_RUN "$INSTDIR\visit.exe"
!insertmacro MUI_PAGE_FINISH

; Uninstaller pages
!insertmacro MUI_UNPAGE_INSTFILES

; Language files
!insertmacro MUI_LANGUAGE "English"


; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "..\installation\visit${PRODUCT_VERSION}.exe"
InstallDir "$PROGRAMFILES\LLNL\VisIt ${PRODUCT_VERSION}"
InstallDirRegKey HKLM "${PRODUCT_DIR_REGKEY}" ""
ShowInstDetails show
ShowUnInstDetails show

Var CreatedPythonLinks
Var DefaultDatabase
Var SelectingDefaultDatabase
Var InstallDevelopmentFiles

###############################################################################
#
# Functions
#
###############################################################################

Function .onInit
  ;Extract InstallOptions INI files
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "NetworkConfig.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "ParallelBank.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "WantDatabasePlugin.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "DefaultDatabasePlugin.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "ChooseInstallDevelopmentFiles.ini"
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "ClickInstall.ini"
  Strcpy $SelectingDefaultDatabase "no"
  Strcpy $DefaultDatabase ""
  Strcpy $InstallDevelopmentFiles "yes"
FunctionEnd

#
# This function is called when we show the Network configuration screen.
#
Function ChooseNetworkConfig
  !insertmacro MUI_HEADER_TEXT "Network configuration" "Select the desired network configuration."
  !insertmacro MUI_INSTALLOPTIONS_DISPLAY "NetworkConfig.ini"
FunctionEnd

#
# This function is called when we show the Choose parallel bank screen.
#
Function ChooseParallelBank
  !insertmacro MUI_HEADER_TEXT "Parallel bank selection" "Select the desired parallel bank."
  !insertmacro MUI_INSTALLOPTIONS_DISPLAY "ParallelBank.ini"
FunctionEnd

#
# This function is called when we show the WantDatabasePlugin screen. We read whether
# the user chose yes or no and use that to set a variable that we use to determine the
# look of the next page.
#
Function WantDefaultDatabasePlugin
  !insertmacro MUI_HEADER_TEXT "Select default database reader plugin" "Do you want to select a default database reader plugin?"
  !insertmacro MUI_INSTALLOPTIONS_DISPLAY "WantDatabasePlugin.ini"

  # Get whether or not the user selected a default file format.
  !insertmacro MUI_INSTALLOPTIONS_READ $0 "WantDatabasePlugin.ini" "Field 1" "State"
  # If $0=="1" then we're going to have a networkconfig
  Strcmp $0 "1" NoDefaultDatabase PickedDefaultDatabase
PickedDefaultDatabase:
    # We got here because we picked a default database. Enable the database plugin combobox.
    Strcpy $SelectingDefaultDatabase "yes"
    Goto EndWantDefaultDatabasePlugin
NoDefaultDatabase:
    # We got here because we picked no default database.
    Strcpy $SelectingDefaultDatabase "no"
EndWantDefaultDatabasePlugin:
FunctionEnd

#
# This function is called when we want to actually choose a database reader plugin. If the
# user chose not to set up a database reader plugin then we show the "Click install"
# screen instead.
#
Function ChooseDefaultDatabasePlugin
  Strcmp $SelectingDefaultDatabase "yes" YesPickDatabase NoPickDatabase
YesPickDatabase:
  !insertmacro MUI_HEADER_TEXT "Select default database reader plugin" "Select the database reader plugin that VisIt will try first when opening a database."
  !insertmacro MUI_INSTALLOPTIONS_DISPLAY "DefaultDatabasePlugin.ini"
  !insertmacro MUI_INSTALLOPTIONS_READ $0 "DefaultDatabasePlugin.ini" "Field 1" "State"
   Strcpy $DefaultDatabase "-default_format $0"
NoPickDatabase:
FunctionEnd

#
# Allow the user to choose whether import libraries and header files get installed.
#
Function ChooseInstallDevelopmentFiles
    !insertmacro MUI_HEADER_TEXT "Install plugin development files" "Do you want to want to install files for plugin development?"
    !insertmacro MUI_INSTALLOPTIONS_DISPLAY "ChooseInstallDevelopmentFiles.ini"
    
    # Get whether or not the user wanted plugin development files
    !insertmacro MUI_INSTALLOPTIONS_READ $0 "ChooseInstallDevelopmentFiles.ini" "Field 1" "State"
    # If $0=="1" then we're going to install plugin development files
    Strcmp $0 "1" NoInstallDevelopmentFiles YesInstallDevelopmentFiles
NoInstallDevelopmentFiles:
    # We got here because we did not want to install the "plugin dev" section
    SectionSetFlags 7 0
    Strcpy $InstallDevelopmentFiles "no"
    Goto EndInstallDevelopmentFiles
YesInstallDevelopmentFiles:
    Strcpy $InstallDevelopmentFiles "yes"
EndInstallDevelopmentFiles:
FunctionEnd

###############################################################################
#
# Sections
#
###############################################################################
Section "Executable Components" SEC01
  SetOutPath "$INSTDIR"
  SetOverwrite ifnewer
  File "..\bin\MSVC7.Net\Release\*.dll"
  File "..\bin\MSVC7.Net\Release\*.exe"
  File "..\bin\MSVC7.Net\Release\visit-config-*.ini"
  File "..\bin\MSVC7.Net\Release\xml2plugin.bat"
  File "..\bin\MSVC7.Net\Release\makemovie.py"
  File "..\bin\MSVC7.Net\Release\makemoviemain.py"

  # Icon files
  File "..\resources\*.ico"
  # Qt DLL
  File "${QTPATH}\lib\${QTDLL}"
  # MSVC 7 .NET runtime libraries
  File "C:\Program files\Microsoft Visual Studio .NET 2003\SDK\v1.1\Bin\msvcr71.dll"
  File "C:\Program files\Microsoft Visual Studio .NET 2003\SDK\v1.1\Bin\msvcp71.dll"
SectionEnd

Section "Database plugins" SEC02
  SetOutPath "$INSTDIR\databases"
  File "..\bin\MSVC7.Net\Release\databases\lib*.dll"
SectionEnd

Section "Plot plugins" SEC03
  SetOutPath "$INSTDIR\plots"
  File "..\bin\MSVC7.Net\Release\plots\lib*.dll"
SectionEnd

Section "Operator plugins" SEC04
  SetOutPath "$INSTDIR\operators"
  File "..\bin\MSVC7.Net\Release\operators\lib*.dll"
SectionEnd

Section "Python modules" SEC05
  SetOutPath "$INSTDIR"
  File /r "..\bin\MSVC7.Net\Release\Python"
SectionEnd

Section HelpFiles SEC06
  SetOutPath "$INSTDIR\help"
  File "..\bin\MSVC7.Net\Release\help\*.html"
  File "..\bin\MSVC7.Net\Release\help\visit.helpml"
  SetOutPath "$INSTDIR"
SectionEnd

Section DataFiles SEC07
  SetOutPath "$INSTDIR\data"
  #
  # This references Files that are on my local C:\ drive since I don't want to have to
  # make the source distribution have projects to build the test data programs.
  #
  File "..\..\VisItData\*.silo"
  File "..\..\VisItData\wave.visit"
  File "..\..\VisItData\PDB\db*.pdb"
  File "..\..\VisItData\ANALYZE_test_data\*.hdr"
  File "..\..\VisItData\ANALYZE_test_data\*.img"
  File "..\..\VisItData\ANALYZE_test_data\*.visit"
  File "..\..\VisItData\FVCOM\*.nc"
  File "..\..\VisItData\molecules\crotamine.pdb"
  File "..\..\VisItData\molecules\1NTS.pdb"
  File "..\..\VisItData\molecules\1UZ9.pdb"
SectionEnd

Section "Plugin development" SEC08
  SetOutPath "$INSTDIR\lib"
  File "..\lib\MSVC7.Net\Release\*.lib"
  SetOutPath "$INSTDIR"
  File /r "..\include"
SectionEnd

Section MyImageDirectory
  # This will hopefully create an image storage directory that VisIt can use.
  SetOutPath "$INSTDIR\My images"

  # Make sure that we're in the VisIt directory by default when the
  # application runs.
  SetOutPath "$INSTDIR"
SectionEnd

Section MovieTemplates
  # Copy over the VisIt movie templates.
  SetOutPath "$INSTDIR\movietemplates"

  File "..\bin\MSVC7.Net\Release\movietemplates\*.mt"
  File "..\bin\MSVC7.Net\Release\movietemplates\*.py"
  File "..\bin\MSVC7.Net\Release\movietemplates\*.ui"
  File "..\bin\MSVC7.Net\Release\movietemplates\*.bmp"
SectionEnd

Section AddVisItRegKeys
#
# This section installs the VISIT<version> key, which tells visit.exe where
# to find the rest of the VisIt components. Note that we put keys in 
# HKEY_LOCAL_MACHINE and in HKEY_CURRENT_USER.
#
  WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "" ""
  WriteRegStr HKCU "VISIT${PRODUCT_VERSION}" "" ""
  WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "VISITHOME" "$INSTDIR"
  WriteRegStr HKCU "VISIT${PRODUCT_VERSION}" "VISITHOME" "$INSTDIR"

  # Write the system config that the user chose.
  !insertmacro MUI_INSTALLOPTIONS_READ $0 "NetworkConfig.ini" "Field 1" "State"
  # If $0=="" then we're going to have a networkconfig
  Strcmp $0 "0" HaveNetworkConfig SkipNetworkConfig
HaveNetworkConfig:
    !insertmacro MUI_INSTALLOPTIONS_READ $0 "NetworkConfig.ini" "Field 2" "State"
    # If $0=="" then we're going to use the closed config
    Strcmp $0 "1" If_OpenNetworkConfig Else_If_ClosedNetworkConfig
If_OpenNetworkConfig:
    Strcpy $r0 "visit-config-open.ini"
    Goto EndIfNetworkConfig
Else_If_ClosedNetworkConfig:
    !insertmacro MUI_INSTALLOPTIONS_READ $0 "NetworkConfig.ini" "Field 3" "State"
    # If $0=="" then we're going to use the closed config
    Strcmp $0 "1" If_ClosedNetworkConfig Else_If_NERSCNetworkConfig
    If_ClosedNetworkConfig:
        Strcpy $r0 "visit-config-closed.ini"
        goto EndIfNetworkConfig
    Else_If_NERSCNetworkConfig:
        !insertmacro MUI_INSTALLOPTIONS_READ $0 "NetworkConfig.ini" "Field 4" "State"
        # If $0=="" then we're going to use the NERSC config
        Strcmp $0 "1" If_NERSCNetworkConfig Else_If_ORNLNetworkConfig
        If_NERSCNetworkConfig:
            Strcpy $r0 "visit-config-nersc.ini"
            goto EndIfNetworkConfig
        Else_If_ORNLNetworkConfig:
            !insertmacro MUI_INSTALLOPTIONS_READ $0 "NetworkConfig.ini" "Field 5" "State"
            # If $0=="" then we're going to use the ORNL config
            Strcmp $0 "1" If_ORNLNetworkConfig EndIfNetworkConfig
            If_ORNLNetworkConfig:
                Strcpy $r0 "visit-config-ornl.ini"
                goto EndIfNetworkConfig
EndIfNetworkConfig:
    # Store the string from $r0 minus the ".ini" extension in $r2
    Strlen $r1 $r0
    IntOp $r1 $r1 - 4
    StrCpy $r2 $r0 $r1
    WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "VISITSYSTEMCONFIG" $r2
    WriteRegStr HKCU "VISIT${PRODUCT_VERSION}" "VISITSYSTEMCONFIG" $r2
    
    # If the user's chosen bank was not bdivp then replace bdivp in the config files
    Strcmp $r2 "bdivp" NetworkConfigDone ReplaceBank
ReplaceBank:
    !insertmacro MUI_INSTALLOPTIONS_READ $r3 "ParallelBank.ini" "Field 1" "State"
    Push $r0
    Push "bdivp"
    Push $r3
    VIKit::ReplaceStringInFile
NetworkConfigDone:
SkipNetworkConfig:

  # Write any additional arguments, like the default database format, to the VISITARGS key.
  Strcmp $SelectingDefaultDatabase "yes" HaveDefaultDatabase NoDefaultDatabase
HaveDefaultDatabase:
     WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "VISITARGS" $DefaultDatabase
     WriteRegStr HKCU "VISIT${PRODUCT_VERSION}" "VISITARGS" $DefaultDatabase
NoDefaultDatabase:

  # If the Python installation path for Python 2.3 does not exist then create it.
#  ReadRegStr $CreatedPythonLinks HKCU "Software\Python\PythonCore\2.3\InstallPath"
#  Strcmp $CreatedPythonLinks "" CreatePythonLinks SkipAddingVisItKeys
#CreatePythonLinks:
  WriteRegStr HKCU "Software\Python\PythonCore\2.3\InstallPath" "" "$INSTDIR"
  WriteRegStr HKCU "Software\Python\PythonCore\2.3\PythonPath"  "" "$INSTDIR\Python\Lib;$INSTDIR\Python\DLLs;$INSTDIR\Python\Lib\lib-tk"
#  WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "CreatedPythonLinks"  "yes"
#SkipAddingVisItKeys:
SectionEnd

Section CreateLinks
  CreateDirectory "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}"
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\VisIt ${PRODUCT_VERSION}.lnk"     "$INSTDIR\visit.exe" ""     "" 0 SW_SHOWMINIMIZED "" "VisIt allows you to visualize simulation data."
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\VisIt ${PRODUCT_VERSION} in stereo.lnk" "$INSTDIR\visit.exe" "-stereo"     "" 0 SW_SHOWMINIMIZED "" "VisIt allows you to visualize simulation data in stereo."
  CreateShortCut "$DESKTOP\VisIt ${PRODUCT_VERSION}.lnk"                                 "$INSTDIR\visit.exe" ""     "" 0 SW_SHOWMINIMIZED "" "VisIt allows you to visualize simulation data."
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\VisIt Command Line Interface.lnk" "$INSTDIR\visit.exe" "-cli" "" 0 SW_SHOWNORMAL    "" "VisIt's command line interface allows you to visualize simulation data via Python scripting."
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Silex.lnk"                        "$INSTDIR\silex.exe" ""     "" 0 SW_SHOWNORMAL    "" "Silex allows you to browse the contents of Silo files."

  # Optionally add a link for xmledit.
  Strcmp $InstallDevelopmentFiles "yes" YesInstallDevelopmentFiles NoInstallDevelopmentFiles
YesInstallDevelopmentFiles:
      CreateDirectory "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\XML Edit.lnk" "$INSTDIR\xmledit.exe" "" "" 0 SW_SHOWNORMAL    "" "XMLEdit allows you to edit the XML files that describe VisIt's plugins."
      CreateDirectory "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation\VTK classes.lnk"    "http://www.vtk.org/doc/release/5.0/html/classes.html"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation\Qt classes.lnk"     "http://doc.trolltech.com/3.0/classes.html"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation\OpenGL library.lnk" "http://www.rush3d.com/reference/opengl-bluebook-1.0/"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation\Python library.lnk" "http://docs.python.org/lib/lib.html"
      CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Plugin development\Documentation\HDF5 library.lnk"   "http://hdf.ncsa.uiuc.edu/HDF5/doc/RM_H5Front.html"
NoInstallDevelopmentFiles:
SectionEnd

Section AddFileAssociations
  # Associate the Silo file format with VisIt and Silex.
  WriteRegStr HKCR ".silo" "" "siloFile"
  WriteRegStr HKCR "siloFile" "" "Silo File"
  WriteRegStr HKCR "siloFile\DefaultIcon" "" "$INSTDIR\silo.ico"
  WriteRegStr HKCR "siloFile\shell\Explore\command" "" '$INSTDIR\silex.exe "%1"'
  WriteRegStr HKCR "siloFile\shell\open\command" "" '$INSTDIR\visit.exe -o "%1"'

  # Associate the VisIt file format with VisIt.
  WriteRegStr HKCR ".visit" "" "visitFile"
  WriteRegStr HKCR "visitFile" "" "VisIt File"
  WriteRegStr HKCR "visitFile\DefaultIcon" "" "$INSTDIR\visitfile.ico"
  WriteRegStr HKCR "visitFile\shell\open\command" "" '$INSTDIR\visit.exe -o "%1"'

  # Associate the VisIt session file format with VisIt.
  WriteRegStr HKCR ".vses" "" "visitSessionFile"
  WriteRegStr HKCR "visitSessionFile" "" "VisIt Session File"
  WriteRegStr HKCR "visitSessionFile\DefaultIcon" "" "$INSTDIR\visitsessionfile.ico"
  WriteRegStr HKCR "visitSessionFile\shell\Make movie\command"          "" '$INSTDIR\visit.exe -movie -format tiff -sessionfile "%1"'
  WriteRegStr HKCR "visitSessionFile\shell\Make 480x480 movie\command"  "" '$INSTDIR\visit.exe -movie -format tiff -geometry 480x480 -sessionfile "%1"'
  WriteRegStr HKCR "visitSessionFile\shell\Make 640x480 movie\command"  "" '$INSTDIR\visit.exe -movie -format tiff -geometry 640x480 -sessionfile "%1"'
  WriteRegStr HKCR "visitSessionFile\shell\Make 800x600 movie\command"  "" '$INSTDIR\visit.exe -movie -format tiff -geometry 800x600 -sessionfile "%1\'
  WriteRegStr HKCR "visitSessionFile\shell\Make 1024x768 movie\command" "" '$INSTDIR\visit.exe -movie -format tiff -geometry 1024x768 -sessionfile "%1"'
  WriteRegStr HKCR "visitSessionFile\shell\Edit\command" "" 'notepad.exe "%1"'
  WriteRegStr HKCR "visitSessionFile\shell\open\command" "" '$INSTDIR\visit.exe -sessionfile "%1"'
SectionEnd

Section AddJavaInstallPath
   # Call our VIkit DLL to get the $INSTDIR variable formatted as a Java preference.
   VIkit::GetInstallPathFormattedForJava
   Pop $R0
   # Write the reformatted string as a Java preference.
   WriteRegStr HKLM "SOFTWARE\JavaSoft\Prefs\llnl\visit" "/V/I/S/I/T/H/O/M/E" $R0
SectionEnd

Section MigrateConfigFiles
   # Call our VIkit DLL to get see if there are config files to migrate.
   Push ${PRODUCT_VERSION}
   VIkit::GetPathToOlderConfigFiles
   # Install path
   Pop $R1
   # Message box prompt
   Pop $R0

   # If $R1 == "" then no config files
   Strcmp $R1 "" NoMoveConfigs HaveConfigs
HaveConfigs:
   # Ask if the user wants to move the configs
   MessageBox MB_YESNO $R0 IDYES YesMoveConfigs IDNO NoMoveConfigs
YesMoveConfigs:
   # Call our VIkit DLL to migrate the config files.
   Push $R1
   VIKit::MigrateConfigFiles
NoMoveConfigs:
SectionEnd

Section -AdditionalIcons
  WriteIniStr "$INSTDIR\${PRODUCT_NAME}.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\VisIt Home Page.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  CreateShortCut "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}\Uninstall VisIt ${PRODUCT_VERSION}.lnk" "$INSTDIR\uninst.exe"
SectionEnd

Section -Post
  WriteUninstaller "$INSTDIR\uninst.exe"
  WriteRegStr HKLM "${PRODUCT_DIR_REGKEY}" "" "$INSTDIR\visit.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayName" "$(^Name)"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString" "$INSTDIR\uninst.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayIcon" "$INSTDIR\visit.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayVersion" "${PRODUCT_VERSION}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "URLInfoAbout" "${PRODUCT_WEB_SITE}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "Publisher" "${PRODUCT_PUBLISHER}"
SectionEnd


Function un.onUninstSuccess
  HideWindow
  MessageBox MB_ICONINFORMATION|MB_OK "$(^Name) was successfully removed from your computer."
FunctionEnd

Function un.onInit
  MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove $(^Name) and all of its components?" IDYES +2
  Abort
FunctionEnd

Section Uninstall
  # Remove the desktop shortcut
  Delete "$DESKTOP\VisIt ${PRODUCT_VERSION}.lnk"

  # Remove the Start menu program group
  RMDir /r "$SMPROGRAMS\VisIt ${PRODUCT_VERSION}"

  # Remove all of the VisIt software components
  RMDir /r "$INSTDIR"

  # Delete the Silo file type from the registry.
  DeleteRegKey HKCR ".silo"
  DeleteRegKey HKCR "siloFile"
  # Delete the VisIt session file type from the registry.
  DeleteRegKey HKCR ".vses"
  DeleteRegKey HKCR "visitSessionFile"
  # Delete the VisIt file type from the registry.
  DeleteRegKey HKCR ".visit"
  DeleteRegKey HKCR "visitFile"

  # If we created links for Python, remove them when we remove VisIt.
#  ReadRegStr $CreatedPythonLinks HKCR "VISIT${PRODUCT_VERSION}" "CreatedPythonLinks"
#  Strcmp $CreatedPythonLinks "yes" RemovePythonLinks SkipRemovingPythonLinks
#RemovePythonLinks:
#  DeleteRegKey HKCU "Software\Python\PythonCore\2.3\InstallPath"
#  DeleteRegKey HKCU "Software\Python\PythonCore\2.3\PythonPath"
#SkipRemovingPythonLinks:

  # Delete the VisIt <version> key registry.
  DeleteRegKey HKCR "VISIT${PRODUCT_VERSION}"
  DeleteRegKey HKCU "VISIT${PRODUCT_VERSION}"

  DeleteRegKey ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}"
  DeleteRegKey HKLM "${PRODUCT_DIR_REGKEY}"

  DeleteRegKey HKLM "SOFTWARE\JavaSoft\Prefs\llnl"

  SetAutoClose true
SectionEnd
