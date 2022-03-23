; -*- compile-command:"/cygdrive/c/Program\ Files/NSIS/makensis.exe xcas.nsi" -*- -*-
;NSIS Modern User Interface version 1.70
;Basic Example Script
;Written by Joost Verburg
; modified for xcas by B. Parisse

;--------------------------------
;Include Modern UI

  !include "MUI.nsh"

;--------------------------------
;General

  ;Name and file
  Name "Xcas"
  OutFile "xcasinst.exe"

  ;Default installation folder
  InstallDir "c:\xcas"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\xcas" ""

;--------------------------------
;Variables

  Var MUI_TEMP
  Var STARTMENU_FOLDER

;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING

;--------------------------------
;Pages
LangString bienvenue_title 1036 "Bienvenue dans le programme d'installation de Xcas"
LangString bienvenue_title 1033 "Welcome to Xcas install program"
; espagnol 
LangString bienvenue_title 1034 "Welcome to Xcas install program"
; grec 
LangString bienvenue_title 1032 "Welcome to Xcas install program"
LangString bienvenue_title 2052 "Welcome to Xcas install program"
!define MUI_WELCOMEPAGE_TITLE $(bienvenue_title)
LangString bienvenue_text 1036 "Cet assistant va installer Xcas, un logiciel libre de calcul formel, de géométrie 2-d et 3-d et tableur. Ce logiciel est publié sous la licence GNU GPL version 3 ou ultérieure, dont vous devrez approuver le contrat juridique rédigé en anglais à la page suivante."
LangString bienvenue_text 1033 "The installer will install Xcas, a free computer algebra system with 2-d/3-d geometry and spreadsheet. Xcas is available under the GNU GPL version 3 or later license"
; espagnol
LangString bienvenue_text 1034 "The installer will install Xcas, a free computer algebra system with 2-d/3-d geometry and spreadsheet. Xcas is available under the GNU GPL version 3 or later license"
; grec
LangString bienvenue_text 1032 "The installer will install Xcas, a free computer algebra system with 2-d/3-d geometry and spreadsheet. Xcas is available under the GNU GPL version 3 or later license"
LangString bienvenue_text 2052 "The installer will install Xcas, a free computer algebra system with 2-d/3-d geometry and spreadsheet. Xcas is available under the GNU GPL version 3 or later license"
!define MUI_WELCOMEPAGE_TEXT $(bienvenue_text)

!insertmacro MUI_PAGE_WELCOME

  !insertmacro MUI_PAGE_LICENSE "c:\xcas\COPYING"
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "xcas" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Menu demarrer"
  
  !insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER

  !define MUI_FINISHPAGE_TITLE "Installation OK"
LangString aurevoir 1036 "Pour lancer Xcas, cliquez sur frxcas dans le menu Démarrer -> Xcas. Vous pouvez aussi cliquer sur le raccourci vers xcasfr.bat sur le bureau.\r\n\r\nVous pouvez aussi utiliser le fichier xcaskey.bat pour lancer xcas s'il est installe sur une clef USB apres avoir modifie le fichier runxcasp.fr qui fonctionne avec un lecteur reseau ou le repertoire personnel est P:.\r\nSi vous etes administrateur, vous pouvez specifier le repertoire home des utilisateurs en vous inspirant de xcasfrjp.bat et runxcasj.fr (par exemple en changeant la lettre du lecteur reseau p dans export XCAS_HOME=/cygdrive/p et export XCAS_AUTOSAVE_FOLDER=/cygdrive/p)\r\nPour compiler une application C++ utilisant giac, lisez le fichier README"
LangString aurevoir 1033 "To run Xcas, select Xcasen from the Xcas submenu of the Start menu. If you want to compile a C++ program using giac, please read the README file from the Xcas installation directory."
; espagnol
LangString aurevoir 1034 "To run Xcas, select Xcasen from the Xcas submenu of the Start menu. If you want to compile a C++ program using giac, please read the README file from the Xcas installation directory."
; grec
LangString aurevoir 1032 "To run Xcas, select Xcasen from the Xcas submenu of the Start menu. If you want to compile a C++ program using giac, please read the README file from the Xcas installation directory."
LangString aurevoir 2052 "To run Xcas, select Xcaszh from the Xcas submenu of the Start menu. If you want to compile a C++ program using giac, please read the README file from the Xcas installation directory."
  !define MUI_FINISHPAGE_TEXT $(aurevoir)
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  

;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"
  !insertmacro MUI_LANGUAGE "French"
  !insertmacro MUI_LANGUAGE "Spanish"
  !insertmacro MUI_LANGUAGE "Greek"
  !insertmacro MUI_LANGUAGE "SimpChinese"



;--------------------------------
;Installer Sections

Section "Xcas" SecDummy

  SetOutPath "$INSTDIR"
  
  ;ADD YOUR OWN FILES HERE...
  File /r c:\xcas\*.*

  Exec '"$INSTDIR\win2unix.exe" $INSTDIR'
  WriteRegStr HKCR ".xws" "" "Xcas.Script"
  WriteRegStr HKCR "Xcas.Script" "" "Xcas Script File"
  WriteRegStr HKCR "Xcas.Script\DefaultIcon" "" "$INSTDIR\xcas.ico"
  WriteRegStr HKCR "Xcas.Script\shell" "" "open"
  WriteRegStr HKCR "Xcas.Script\shell\open\command" "" '$INSTDIR\xcasfr.bat "%1"'

  ;Store installation folder
  WriteRegStr HKCU "xcas" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
  ;Create shortcuts
  CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\update.lnk" "$INSTDIR\update.bat"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\xcasEn.lnk" "$INSTDIR\xcasen.bat"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\xcasGre.lnk" "$INSTDIR\xcasgre.bat"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\xcasEs.lnk" "$INSTDIR\xcases.bat"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\xcasZh.lnk" "$INSTDIR\xcaszh.bat"
  CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\frXcas.lnk" "$INSTDIR\xcasfr.bat"
  ; CreateDirectory "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas" 
  ; CreateShortCut "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  ; CreateShortCut "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\Xcasfr.lnk" "$INSTDIR\xcasfr.bat"
  ; CreateShortCut "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\XcasGre.lnk" "$INSTDIR\xcasgre.bat"
  CreateShortCut "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\Xcasdoc.lnk" "$INSTDIR\doc\index.html"  

  !insertmacro MUI_STARTMENU_WRITE_END
StrCmp $LANGUAGE ${LANG_FRENCH} 0 +2
  CreateShortCut "$DESKTOP\xcasfr.lnk" "$INSTDIR\xcasfr.bat"
StrCmp $LANGUAGE ${LANG_GREEK} 0 +2
  CreateShortCut "$DESKTOP\xcasGre.lnk" "$INSTDIR\xcasgre.bat"
StrCmp $LANGUAGE ${LANG_SIMPCHINESE} 0 +2
  CreateShortCut "$DESKTOP\xcasZh.lnk" "$INSTDIR\xcaszh.bat"
  CreateShortCut "$DESKTOP\xcas_doc.lnk" "$INSTDIR\doc\index.html"

	; Delete previous index files if they exist
  Delete "$INSTDIR\doc\fr\casflan\html_mtt"
  Delete "$INSTDIR\doc\fr\casflan\html_mall"
  Delete "$INSTDIR\doc\fr\casflan\html_vall"
  Delete "$INSTDIR\doc\fr\html_mtt"
  Delete "$INSTDIR\doc\fr\html_mall"
  Delete "$INSTDIR\doc\fr\html_vall"

;  IfFileExists "c:\cygwin" +4 0
;    CreateDirectory "c:\cygwin"
;    CreateDirectory "c:\cygwin\bin"
;    CopyFiles "$INSTDIR\sh.exe" "c:\cygwin\bin"

SectionEnd

;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_SecDummy ${LANG_ENGLISH} "Xcas."
  LangString DESC_SecDummy ${LANG_FRENCH} "Xcas."
  LangString DESC_SecDummy ${LANG_SPANISH} "Xcas."
  LangString DESC_SecDummy ${LANG_SIMPCHINESE} "Xcas."

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${SecDummy} $(DESC_SecDummy)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END

;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;ADD YOUR OWN FILES HERE...

  Delete "$INSTDIR\doc\fr\casflan\*.*"
  RMDir "$INSTDIR\doc\fr\casflan"
  Delete "$INSTDIR\doc\fr\cassim\*.*"
  RMDir "$INSTDIR\doc\fr\cassim"
  Delete "$INSTDIR\doc\fr\castor\*.*"
  RMDir "$INSTDIR\doc\fr\castor"
  Delete "$INSTDIR\doc\fr\casrouge\*.*"
  RMDir "$INSTDIR\doc\fr\casrouge"
  Delete "$INSTDIR\doc\fr\casgeo\*.*"
  RMDir "$INSTDIR\doc\fr\casgeo"
  Delete "$INSTDIR\doc\fr\tutoriel\*.*"
  RMDir "$INSTDIR\doc\fr\tutoriel"
  Delete "$INSTDIR\doc\fr\*.*"
  RMDir "$INSTDIR\doc\fr"
  Delete "$INSTDIR\doc\en\tutoriel_en\*.*"
  RMDir "$INSTDIR\doc\en\tutoriel_en"
  Delete "$INSTDIR\doc\en\*.*"
  RMDir "$INSTDIR\doc\en"
  Delete "$INSTDIR\doc\es\*.*"
  RMDir "$INSTDIR\doc\es"
  Delete "$INSTDIR\doc\*.*"
  RMDir "$INSTDIR\doc\"
  Delete "$INSTDIR\examples\geo\*.*"
  RMDir "$INSTDIR\examples\geo"
  Delete "$INSTDIR\examples\morley\*.*"
  RMDir "$INSTDIR\examples\morley"
  Delete "$INSTDIR\examples\lewisw\*.*"
  RMDir "$INSTDIR\examples\lewisw"
  Delete "$INSTDIR\examples\codage\*.*"
  RMDir "$INSTDIR\examples\codage\"
  Delete "$INSTDIR\examples\polyfact\*.*"
  RMDir "$INSTDIR\examples\polyfact\"
  Delete "$INSTDIR\examples\simulation\*.*"
  RMDir "$INSTDIR\examples\simulation\"
  Delete "$INSTDIR\examples\arit\*.*"
  RMDir "$INSTDIR\examples\arit"
  Delete "$INSTDIR\examples\*.*"
  RMDir "$INSTDIR\examples\"
  Delete "$INSTDIR\locale\es\LC_MESSAGES\*.*"
  RMDir "$INSTDIR\locale\es\LC_MESSAGES"
  RMDir "$INSTDIR\locale\es\"
  Delete "$INSTDIR\locale\de\LC_MESSAGES\*.*"
  RMDir "$INSTDIR\locale\de\LC_MESSAGES"
  RMDir "$INSTDIR\locale\de\"
  Delete "$INSTDIR\locale\el\LC_MESSAGES\*.*"
  RMDir "$INSTDIR\locale\el\LC_MESSAGES"
  RMDir "$INSTDIR\locale\el\"
  Delete "$INSTDIR\locale\fr\LC_MESSAGES\*.*"
  RMDir "$INSTDIR\locale\fr\LC_MESSAGES"
  RMDir "$INSTDIR\locale\fr\"
  Delete "$INSTDIR\locale\*.*"
  RMDir "$INSTDIR\locale\"
  Delete "$INSTDIR\AsTeX\*.*"
  RMDir "$INSTDIR\AsTeX"
  CreateDirectory "$INSTDIR_"
  Rename "$INSTDIR\*.xws" "$INSTDIR_data"
  Rename "$INSTDIR\*.cas" "$INSTDIR_data"
  Rename "$INSTDIR\*.tex" "$INSTDIR_data"
  Rename "$INSTDIR\*.eps" "$INSTDIR_data"
  Rename "$INSTDIR\*.png" "$INSTDIR_data"
  Delete "$INSTDIR\*.*"
  RMDir "$INSTDIR"

  !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    
  Delete "$SMPROGRAMS\$MUI_TEMP\Uninstall.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\update.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\xcasEn.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\xcasEs.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\frXcas.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\XcasGre.lnk"
  Delete "$DESKTOP\xcasfr.lnk"
  ; Delete "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\Xcasfr.lnk"
  ; Delete "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\Uninstall.lnk"
  ; RMDir "C:\Documents and Settings\All Users\Menu Démarrer\Programmes\Xcas\"

  ;Delete empty start menu parent diretories
  StrCpy $MUI_TEMP "$SMPROGRAMS\$MUI_TEMP"
 
  startMenuDeleteLoop:
    RMDir $MUI_TEMP
    GetFullPathName $MUI_TEMP "$MUI_TEMP\.."
    
    IfErrors startMenuDeleteLoopDone
  
    StrCmp $MUI_TEMP $SMPROGRAMS startMenuDeleteLoopDone startMenuDeleteLoop
  startMenuDeleteLoopDone:



  DeleteRegKey /ifempty HKCU "Software\Xcas"

SectionEnd


