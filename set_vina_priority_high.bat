@echo off
REM Sets the running vina.exe process to High CPU priority (WMIC setpriority
REM code 128). Use this to speed up docking at the cost of leaving less CPU
REM for other programs while it runs.
wmic process where name="vina.exe" CALL setpriority 128
