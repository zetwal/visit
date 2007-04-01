; Script generated by the HM NIS Edit Script Wizard.

##############################################################################
#
# File: sourceinstallation.nsi
#
# Purpose: This file contains the instructions that the NSIS installer needs
#          in order to create an installation program for VisIt's source code.
#
# Programmer: Brad Whitlock
# Date:       Mon, Dec 15 11:26:34 PDT 2003
#
# Modifications:
#   Brad Whitlock, Mon Feb 9 14:45:04 PST 2004
#   Updated for 1.2.7
#
##############################################################################

; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "VisIt for Windows source code"
!define PRODUCT_VERSION "1.2.7"
!define PRODUCT_PUBLISHER "LLNL"
!define PRODUCT_WEB_SITE "http://www.llnl.gov/visit"
!define PRODUCT_DIR_REGKEY "Software\Microsoft\Windows\CurrentVersion\App Paths\visitdev${PRODUCT_VERSION}"

SetCompressor bzip2

; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "..\resources\visit.ico"

; Welcome page
!insertmacro MUI_PAGE_WELCOME
; License page
!insertmacro MUI_PAGE_LICENSE "copyright.txt"
; Directory page
!insertmacro MUI_PAGE_DIRECTORY
; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\VisItBuildInstructionsOnWindows.doc"
!insertmacro MUI_PAGE_FINISH

; Language files
!insertmacro MUI_LANGUAGE "English"

; Reserve files
!insertmacro MUI_RESERVEFILE_INSTALLOPTIONS

; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "..\installation\visitdev${PRODUCT_VERSION}.exe"
InstallDir "C:\VisItDev${PRODUCT_VERSION}"
InstallDirRegKey HKLM "${PRODUCT_DIR_REGKEY}" ""
ShowInstDetails show

Section "ProjectFiles" SEC02
  SetOutPath "$INSTDIR"
  File /r "..\projects"
SectionEnd

Section "IncludeFiles" SEC03
  SetOutPath "$INSTDIR"
  File /r "..\include"
SectionEnd

Section "InstallationFiles" SEC04
  # Get the instructions file
  SetOutPath "$INSTDIR"
  File "..\VisItBuildInstructionsOnWindows.doc"
  File "..\BUILD_NOTES.txt"
  # The installation files that we use to build the binary and source distributions.
  SetOutPath "$INSTDIR\installation"
  File "binaryinstallation.nsi"
  File "sourceinstallation.nsi"
  File "copyright.txt"
  File "NetworkConfig.ini"
SectionEnd

Section "ScriptFiles" SEC05
  # The script files that we need to build on Windows.
  SetOutPath "$INSTDIR"
  File /r "..\script"
SectionEnd

Section "BinFiles" SEC06
  SetOutPath "$INSTDIR"
  File /r "..\bin"
SectionEnd

Section "VisItSource" SEC07
  SetOutPath "$INSTDIR"
  File /r "..\visit"
SectionEnd

Section "LibFiles" SEC08
  SetOutPath "$INSTDIR"
  File /r "..\lib"
SectionEnd

Section "Resources" SEC09
  SetOutPath "$INSTDIR"
  File /r "..\resources"
SectionEnd

Section -Post
  WriteRegStr HKLM "${PRODUCT_DIR_REGKEY}" "" "$INSTDIR\makensis.exe"
  
  # Set the VISITDEVDIR key in the registry.
  WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "VISITDEVDIR" "$INSTDIR"

  # Set the VISITDEVDIR key in the registry so it will be set as an 
  # environment variable for the current user.
  WriteRegStr HKCU "Environment" "VISITDEVDIR" "$INSTDIR"
SectionEnd
