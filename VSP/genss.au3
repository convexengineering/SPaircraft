#cs ----------------------------------------------------------------------------

 AutoIt Version: 3.3.14.2
 Author:         myName

 Script Function:
	Template AutoIt script.

#ce ----------------------------------------------------------------------------

; Script Start - Add your code below here

#include <MsgBoxConstants.au3>
#include <AutoItConstants.au3>
MsgBox($MB_OK, "Tutorial", "Are you sure you are ready for this?!")

MouseMove(-1520,-20,20)

MouseClick($MOUSE_CLICK_LEFT)

MouseMove(-1520,170,20)

MouseClick($MOUSE_CLICK_LEFT)

Send("reload.vspscript")

MouseMove(-1370, 740, 100)

MouseClick($MOUSE_CLICK_LEFT)

; Taking screenshot

MouseMove(-1420,-20,100)

MouseClick($MOUSE_CLICK_LEFT)

MouseMove(-1420,100,100)

MouseClick($MOUSE_CLICK_LEFT)

MouseMove(-1150,700,20)

MouseClick($MOUSE_CLICK_LEFT)

Send("ss1.png")

MouseMove(-1370, 740, 20)

MouseClick($MOUSE_CLICK_LEFT)

Send("!{F4}")

