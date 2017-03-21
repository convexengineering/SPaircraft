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

For $i = 0 to 8 Step 1

	; Changing the value within reload.vspscript
	MouseMove(-446,121)
	MouseClick($MOUSE_CLICK_LEFT)
	Send("{BACKSPACE}")
	Send($i)
	Send("^s")

	; Getting to Run Script window
	MouseMove(-1520,-20,20)
	MouseClick($MOUSE_CLICK_LEFT)
	MouseMove(-1520,170,20)
	MouseClick($MOUSE_CLICK_LEFT)

	; Clicking reload.vspscript
	MouseMove(-1500,135,20)
	MouseClick($MOUSE_CLICK_LEFT)
	MouseMove(-1370, 740, 20)
	MouseClick($MOUSE_CLICK_LEFT)
	Sleep(1500)

	; Getting right orientation and zoom
	Send("{F5}")
	Send("f")

	; Taking screenshot
	MouseMove(-1420,-20,20)
	MouseClick($MOUSE_CLICK_LEFT)
	MouseMove(-1420,100,20)
	MouseClick($MOUSE_CLICK_LEFT)
	MouseMove(-1150,700,20)
	MouseClick($MOUSE_CLICK_LEFT)
	MouseMove(-1360,650,20)
	MouseClick($MOUSE_CLICK_LEFT)
	Send($i)
	Send(".png")
	; Accept screenshot
	MouseMove(-1370, 740, 20)
	MouseClick($MOUSE_CLICK_LEFT)

	Send("!{F4}")
Next

