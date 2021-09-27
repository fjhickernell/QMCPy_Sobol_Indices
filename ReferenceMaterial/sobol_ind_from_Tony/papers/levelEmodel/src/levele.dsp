# Microsoft Developer Studio Project File - Name="levele" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=levele - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "levele.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "levele.mak" CFG="levele - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "levele - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "levele - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "levele - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "levele - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "levele - Win32 Release"
# Name "levele - Win32 Debug"
# Begin Source File

SOURCE=.\batman.f
# End Source File
# Begin Source File

SOURCE=.\bios.f
DEP_F90_BIOS_=\
	".\cbios.h"\
	".\cnucld.h"\
	".\coutp.h"\
	".\cparam.h"\
	".\cspop.h"\
	".\ctime.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\coeff.f
# End Source File
# Begin Source File

SOURCE=.\gtm1.f
DEP_F90_GTM1_=\
	".\cbios.h"\
	".\cfarf.h"\
	".\cnucld.h"\
	".\coutp.h"\
	".\cparam.h"\
	".\ctime.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\levele.f
DEP_F90_LEVEL=\
	".\cfast.h"\
	".\cparam.h"\
	".\cspop.h"\
	".\ctime.h"\
	".\cvaria.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\levele_model.f
DEP_F90_LEVELE=\
	".\cbios.h"\
	".\cfarf.h"\
	".\cname.h"\
	".\cnearf.h"\
	".\cnucld.h"\
	".\const.h"\
	".\coutp.h"\
	".\cparam.h"\
	".\cspop.h"\
	".\ctime.h"\
	".\cvaria.h"\
	".\par.h"\
	".\rel.h"\
	".\var.h"\
	
# End Source File
# Begin Source File

SOURCE=.\mgspac.f
DEP_F90_MGSPA=\
	".\cfarf.h"\
	".\cparam.h"\
	".\cvaria.h"\
	".\par.h"\
	".\var.h"\
	
# End Source File
# Begin Source File

SOURCE=.\mgtime.f
DEP_F90_MGTIM=\
	".\cfarf.h"\
	".\cnearf.h"\
	".\cnucld.h"\
	".\cparam.h"\
	".\ctime.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\nearf.f
DEP_F90_NEARF=\
	".\cnearf.h"\
	".\cnucld.h"\
	".\coutp.h"\
	".\cparam.h"\
	".\ctime.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\setzer.f
DEP_F90_SETZE=\
	".\cnucld.h"\
	".\coutp.h"\
	".\cparam.h"\
	".\cspop.h"\
	".\ctime.h"\
	".\par.h"\
	
# End Source File
# Begin Source File

SOURCE=.\sift.f
# End Source File
# Begin Source File

SOURCE=.\tres.f
# End Source File
# End Target
# End Project
