@echo off
REM Sets the running vina.exe process to Idle CPU priority (WMIC setpriority
REM code 64). Use this to let docking run in the background without slowing
REM down other programs on the machine.
wmic process where name="vina.exe" CALL setpriority 64
