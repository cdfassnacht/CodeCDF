;+
; NAME:
;       acsproc_process__define
;
; PURPOSE:
;
;       An object for creating a window for a general acsproc process 
;       containing a message window, progress bar, and cancel button.
;     
;
; AUTHORS:
; 
;       showprogress__define was written by:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;       
;       Refactored into a subprogram of acsproc for GUI management
;       By: Joel Brownstein, October, 2010
;
; KEYWORDS:
;
;      CANCELBUTTON: Set this keyword if a Cancel button is desired.
;      DELAY: The total time the widget should be on the display in AutoUpDate
;           mode. The keyword applies only to AutoUpDate mode. Default is 5 seconds.
;      STEPS: The number of steps to take in AutoUpDate mode. The keyword applies only
;           to AutoUpDate mode.
;      message_info: The text of the message_label above the progress bar. Default is "Operation
;           in Progress...".
;      TITLE: ; The text of the top-level base title bar. Default is ""
;      COLOR: The color to draw the progress bar.
;      XSIZE: The XSize of the progress bar in Device coordinates. Default is 150.
;      YSIZE: The YSize of the progress bar in Device coordinates. Default is 10.
;      AUTOUPDATE: Set this keyword to be in AutoUpDate mode.

;
; PROCEDURE:
;       There are two modes. In AutoUpDate mode, a delay and number of steps is
;       required. The modal widget stays on the display until the total time
;       exceeds the DELAY or the requested number of steps is taken. A TIMER
;       widget is used to generate update events. Nothing can be going on
;       concurrently in AutoUpDate mode. To enter AutoUpDate mode, type this:
;
;          progressBar = Obj_New("acsproc_process", /AutoUpDate, Delay=2, Steps=10)
;          progressBar->Start
;          Obj_Destroy, progressBar
;
;       The program will update and destroy itself automatically. (The object
;       itself is not destroyed. You must do this explicitly, as in the example
;       above.)
;
;       In normal mode, the user is responsible for starting, updating, and
;       destroying the progress indicator. The sequence of commands might look
;       like this:
;
;          progressBar = Obj_New("acsproc_process")
;          progressBar->Start
;          FOR j=0,9 DO BEGIN
;             Wait, 0.5  ; Would probably be doing something ELSE here!
;             progressBar->Update, (j+1)*10
;          ENDFOR
;          progressBar->Destroy
;          Obj_Destroy, progressBar
;
;       Normal mode gives you the opportunity to update the Progress Bar
;       in a loop while something else is going on. See the example program
;       at the end of this file.
;
;       Note that the object itself is not destroyed when calling the DESTROY
;       method. You must explicitly destroy the object, as in the example above.
;
; METHODS:
;
;       CHECKCANCEL: This function method returns a 1 if the user has clicked
;           the CANCEL button. Otherwise, it returns a 0.
;
;          cancelled = progressBar->CheckCancel()
;          IF cancelled THEN progressBar->Destroy
;
;       DESTROY: Destroys the acsproc_process widgets. It does NOT destroy the object.
;
;          progressBar->Destroy
;
;       GETPROPERTY: Gets the properties that can be set in the INIT method, including
;          the parent widget ID.
;
;          progressBar->GetProperty, Steps=currentNSteps, Delay=currentDelay
;
;       SETCOLOR: Changes the color of the progress bar.
;
;           progressBar->SetColor, !P.Color
;
;       SETmessage: Changes the text on the widget message.
;
;           progressBar->Setmessage, 'This message instead'
;           progressBar->Setmessage, /clear
;           progressBar->Setmessage, 'Message', /append
;           progressBar->Setmessage, 'Append this after Message', /append
;           progressBar->Setmessage, 'This message label instead', /label
;
;       SETPROPERTY: Allows the user to set the INIT parameter via keywords.
;
;          progressBar->SetProperty, Color=244, XSize=200, message_info='Please Wait...'
;
;       START: Puts the acsproc_process bar on the display. In AutoUpDate mode, the
;          widget starts to automatically update.
;
;          progressBar->Start
;
;       UPDATE: Updates the progress bar. Requires on argument, a number between 0
;          and 100 that indicates the percent of progress bar that should be filled
;          with a color.
;
;          progressBar->Update, 50
;
; EXAMPLE:
;
;       See the example program at http://www.dfanning.com/programs/acsproc_process__define.pro
;
; RESTRICTIONS:
;
;       In contradiction to the IDL documentation, making the parent widget
;          insensitive in normal mode does NOT prevent the parent widgets from
;          receiving events on my Windows NT 4.0, SP 4 system running IDL 5.2,
;          IDL 5.2.1, or IDL 5.3 (beta).
;
;       Note that if you specify a CANCEL button the Show Progress program CANNOT
;       run as a MODAL widget program. Thus, user *may* be able to generate events
;       in the calling program while this program is running.
;
; MODIFICATION HISTORY:
;       Written by:  David Fanning, 26 July 1999.
;       Added code so that the current graphics window doesn't change. 1 September 1999. DWF.
;       Added yet more code for the same purpose. 3 September 1999. DWF.
;       Added a CANCEL button and made other minor modifications. 12 Oct 1999. DWF.
;       
;       Added compartments for process info, progressbar and a buttonbar. 11 Oct 2010. Joel Brownstein
;-


FUNCTION acsproc_process::CheckCancel

; This method checks for a CANCEL Button event. It returns 1
; if an event has occurred and 0 otherwise.

RETURN, self.cancel or not Widget_info(self.tlb, /visible) 

END; -----------------------------------------------------------------------------



PRO acsproc_process::Setmessage, newmessage, label=label, append=append

; This method allows the widget message_label to be changed while
; the program is on the display.
IF keyword_set(label) then begin
  IF N_Elements(newmessage) EQ 0 THEN self.message_label = "" else self.message_label=newmessage
  Widget_Control, self.message_labelID, Set_Value=self.message_label
ENDIF ELSE begin
  IF N_Elements(newmessage) EQ 0 THEN self.message_info = "" else self.message_info=newmessage
  if keyword_set(append) then begin
    Widget_Control, self.message_infoID, get_Value=message_info
    if (n_elements(message_info) eq 1) then if (message_info[0] eq '') then append=0L
  endif
  Widget_Control, self.message_infoID, Set_Value=self.message_info, append=append
ENDELSE

END; -----------------------------------------------------------------------------



PRO acsproc_process::SetCancel, value

; This method sets the Cancel flag.

IF N_Elements(value) EQ 0 THEN value = 1
self.cancel = value
END; -----------------------------------------------------------------------------



PRO acsproc_process::SetColor, color

; This method sets the Cancel flag.

IF N_Elements(color) EQ 0 THEN color = !P.Color
self.color = color
END; -----------------------------------------------------------------------------



PRO acsproc_process::Destroy

; This method takes the widget off the display.

Widget_Control, self.tlb, Destroy=1
self.cancel=0L
END; -----------------------------------------------------------------------------



PRO acsproc_process::UpDate, percent


; This method updates the display. It should be called with
; manual operation. PERCENT should be a value between 0 and 100.

Catch, theError
IF theError NE 0 THEN BEGIN
    Catch, /Cancel

      ; Catch a WSET error silently.

    IF !Error_State.Code EQ -386 THEN RETURN
    message, !Error_State.Msg, /Informational
    RETURN
ENDIF

self.percent = 0 > percent < 100

   ; Update the progress box.

thisWindow = !D.Window
WSet, self.wid
x1 = 0
y1 = 0
x2 = Fix(self.xsize  * (self.percent/100.0))
y2 = self.ysize
Polyfill, [x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], /Device, Color=self.color
IF thisWindow GE 0 AND thisWindow NE self.wid THEN WSet, thisWindow

   ; Check for a CANCEL button event.

IF Widget_Info(self.cancelID, /Valid_ID) THEN BEGIN
   event = Widget_Event(self.cancelID, /NoWait)
   name = Tag_Names(event, /Structure_Name)
   IF name EQ 'WIDGET_BUTTON' THEN self.cancel = 1
ENDIF

END; -----------------------------------------------------------------------------

PRO acsproc_process::Step, step

self->Update, step * 100.0/self.nsteps


END; -----------------------------------------------------------------------------

PRO acsproc_process::Increment, increment

self->Update, self.percent + increment * 100/self.nsteps


END; -----------------------------------------------------------------------------

PRO acsproc_process::Timer_Events, event

; This method processes the widget TIMER events.

Catch, theError
IF theError NE 0 THEN BEGIN
    Catch, /Cancel

      ; Catch a WSET error silently.

    IF !Error_State.Code EQ -386 THEN RETURN
    message_info, !Error_State.Msg, /Informational
    RETURN
ENDIF

   ; Do it the specified number of times, then quit.

self.count = self.count + 1
IF self.count GE self.nsteps THEN BEGIN
   thisWindow = !D.Window
   selfWindow = self.wid
   WSet, self.wid
   Polyfill, [0, 0, self.xsize, self.xsize, 0], $
             [0, self.ysize, self.ysize, 0, 0], /Normal, Color=self.color
   Widget_Control, self.tlb, Destroy=1
   IF thisWindow GE 0 AND thisWindow NE selfWindow THEN WSet, thisWindow
   RETURN
ENDIF

   ; If time gets away from you, then quit.

theTime = Systime(1) - self.startTime
IF theTime GT self.delay THEN BEGIN
   thisWindow = !D.Window
   WSet, self.wid
   Polyfill, [0, 0, self.xsize, self.xsize, 0], $
             [0, self.ysize, self.ysize, 0, 0], /Normal, Color=self.color
   Widget_Control, self.tlb, Destroy=1
   IF thisWindow GE 0 AND thisWindow NE selfWindow THEN WSet, thisWindow
   RETURN
ENDIF

   ; Update the progress box.

thisSize = (self.xsize / self.nsteps) * self.count
thisSize = Float(thisSize) / self.xsize * 100
self->Update, thisSize

   ; Set the next timer event if the CANCEL button is not set.

IF self.cancel EQ 1 THEN BEGIN
   self->Destroy
ENDIF ELSE BEGIN
   Widget_Control, self.tlb, Timer=self.step
ENDELSE
END; -----------------------------------------------------------------------------



PRO acsproc_process_Event, event

; This is the event handler for the program. It simply calls
; the event handling method.

Widget_Control, event.top, Get_UValue=self
thisEvent = Tag_Names(event, /Structure_Name)
IF thisEvent EQ 'WIDGET_BUTTON' THEN self->SetCancel, 1 ELSE $
   self->Timer_Events, event
END; -----------------------------------------------------------------------------


PRO acsproc_process::Start, steps=steps

; This is the method that puts the timer on the display and gets
; things going. The initial timer event is set here.

   ; Find the window index number of any open display window.

thisWindow = !D.Window

   ; Realize the widget.

Widget_Control, self.tlb, /Realize

   ; Set an initial start time.

self.startTime = Systime(1)

   ; Get the window index number of the draw widget.

Widget_Control, self.drawID, Get_Value=wid
self.wid = wid

   ; Back to the open display window.

IF thisWindow GE 0 THEN WSet, thisWindow

if keyword_set(steps) then self.nsteps = steps else self.nsteps = 100

IF self.autoupdate THEN BEGIN

      ; Set the first timer event.

   Widget_Control, self.tlb, Timer=self.step

      ; Register with XMANAGER so you can receive events.

   XManager, 'acsproc_process', self.tlb, Cleanup='acsproc_process_CleanUp'

ENDIF ELSE BEGIN

   IF Widget_Info(self.parent, /Valid_ID) THEN $
      Widget_Control, self.parent, Sensitive=0
   self->Update, 0
   Widget_Control, self.tlb, Kill_Notify='acsproc_process_CleanUp'

ENDELSE

END; -----------------------------------------------------------------------------



PRO acsproc_process_Cleanup, tlb

; This is the cleanup method for the widget. The idea
; here is to reinitialize the widget after the old
; widget has been destroyed. This is necessary, because
; it is not possible to MAP and UNMAP a modal widget.

Widget_Control, tlb, Get_UValue=self
self->ReInitialize

END; -----------------------------------------------------------------------------



PRO acsproc_process::ReInitialize

; This method just reinitializes the acsproc_process widget after the
; previous one has been destroyed. This is called from the widget's
; CLEANUP routine.

IF NOT Float(self.autoupdate) THEN BEGIN
   IF Widget_Info(self.parent, /Valid_ID) THEN $
      Widget_Control, self.parent, Sensitive=1, /Clear_Events
ENDIF

   ; Create the widgets.

IF self.parent EQ -1 THEN BEGIN
   self.tlb = Widget_Base(Title=self.title, Column=1, Base_Align_Center=1)
ENDIF ELSE BEGIN
   IF self.cancelButton THEN modal = 0 ELSE modal = 1
   self.tlb = Widget_Base(Group_Leader=self.parent, Modal=modal, Title=self.title, $
      Column=1, Base_Align_Center=1, Floating=1)
ENDELSE
self.info_baseID = Widget_Base(self.tlb, Column=1, Base_Align_Center=1)
self.acsproc_process_baseID = Widget_Base(self.tlb, Column=3)
self.message_labelID = Widget_label(self.info_baseID, Value=self.message_label, /Dynamic_Resize, /align_left,ysize=24,xsize=750,frame=1)
self.message_infoID = Widget_text(self.info_baseID, Value=self.message_info, ysize=self.message_ysize,xsize=self.message_xsize,/align_left,/scroll)
self.drawID = Widget_Draw(self.acsproc_process_baseID, XSize=self.xsize, YSize=self.ysize,/align_left)
Widget_Control, self.tlb, Set_UValue=self

IF self.cancelButton THEN BEGIN
   self.cancelID = Widget_Label(self.acsproc_process_baseID, Value='                                        ',/align_right)
   self.cancelID = Widget_Button(self.acsproc_process_baseID, Value='Cancel',/align_right)
ENDIF ELSE self.cancelID = -1L

   ; Center the top-level base.

Device, Get_Screen_Size=screenSize
xCenter = screenSize(0) / 2
yCenter = screenSize(1) / 2

geom = Widget_Info(self.tlb, /Geometry)
xHalfSize = geom.Scr_XSize / 2
yHalfSize = geom.Scr_YSize / 2

Widget_Control, self.tlb, XOffset = xCenter-xHalfSize, $
   YOffset = yCenter-yHalfSize

   ; Reset the counter.

self.count = 0
END; -----------------------------------------------------------------------------



PRO acsproc_process::Cleanup

; This CLEANUP method is not usually required, since other widget
; programs are destroying the widget, but is here for completeness.

IF Widget_Info(self.tlb, /Valid_ID) THEN Widget_Control, self.tlb, /Destroy
END; -----------------------------------------------------------------------------



PRO acsproc_process::GetProperty, Parent=parent, Delay=delay, Steps=nsteps, $
   message_info=message_info, Title=title, Color=color, XSize=xsize, $
   YSize=ysize, AutoUpdate=autoupdate

   ; This method allows you to get all the properties available in the INIT method.

parent = self.parent
delay = self.delay
nsteps = self.nsteps
message_label = self.message_label
message_info = self.message_info
message_xsize = self.message_xsize
message_ysize = self.message_ysize
title = self.title
color = self.color
xsize = self.xsize
ysize = self.ysize
autoupdate = self.autoupdate
END; -----------------------------------------------------------------------------



PRO acsproc_process::SetProperty, Parent=parent, Delay=delay, Steps=nsteps, $
   message_label=message_label, message_info=message_info, message_xsize=message_xsize, message_ysize=message_ysize, Title=title, Color=color, XSize=xsize, $
   YSize=ysize, AutoUpdate=autoupdate

   ; This method allows you to set all the properties available in the INIT method.

IF N_Elements(parent) NE 0 THEN BEGIN
    self.parent = parent
   Widget_Control, self.tlb, /Destroy
ENDIF
IF N_Elements(delay) NE 0 THEN self.delay = delay
IF N_Elements(nsteps) NE 0 THEN BEGIN
   self.nsteps = nsteps
   self.step = (Float(self.delay) / self.nsteps)
ENDIF
IF N_Elements(message_label) NE 0 THEN self.message_label = message_label
IF N_Elements(message_info) NE 0 THEN self.message_info = message_info
IF N_Elements(message_xsize) NE 0 THEN self.message_info = message_xsize
IF N_Elements(message_ysize) NE 0 THEN self.message_info = message_ysize
IF N_Elements(title) NE 0 THEN self.title = title
IF N_Elements(color) NE 0 THEN self.color = color
IF N_Elements(xsize) NE 0 THEN self.xsize = xsize
IF N_Elements(ysize) NE 0 THEN self.ysize = ysize
IF N_Elements(autoupdate) NE 0 THEN self.autoupdate = autoupdate
Widget_Control, self.tlb, /Destroy
END; -----------------------------------------------------------------------------



FUNCTION acsproc_process::Init, $
   parent, $             ; The widget ID of the group leader.
   CancelButton=cancelButton, $ ; This keyword is set if a cancel button is required.
   Delay=delay, $        ; The total time the widget should be on the display in AutoUpDate mode.
   Steps=nsteps, $       ; The number of steps to take in AutoUpDate mode.
   message_label=message_label, $    ; The text of the message_label above the progress bar.
   message_info=message_info, $    ; The text of the message_label above the progress bar.
   message_xsize=message_xsize, $    ; The text of the message_label above the progress bar.
   message_ysize=message_ysize, $    ; The text of the message_label above the progress bar.
   Title=title, $        ; The text of the top-level base title bar.
   Color=color, $        ; The color to draw the progress bar.
   XSize=xsize, $        ; The XSize of the progress bar in Device coordinates.
   YSize=ysize, $        ; The YSize of the progress bar in Device coordinates.
   AutoUpdate=autoupdate ; Set this keyword to be in AutoUpDate mode.

   ; A group leader widget (i.e., a parent parameter) is REQUIRED for MODAL operation.

   ; Check keywords.

IF N_Elements(delay) EQ 0 THEN delay = 5
IF N_Elements(nsteps) EQ 0 THEN nsteps = 10
   theStep = (Float(delay) / nsteps)
IF N_Elements(message_label) EQ 0 THEN message_label = '>'
IF N_Elements(message_info) EQ 0 THEN message_info = 'Process Status:'
IF N_Elements(message_xsize) EQ 0 THEN message_xsize = 30
IF N_Elements(message_ysize) EQ 0 THEN message_ysize = 12
IF N_Elements(title) EQ 0 THEN title = ""
IF N_Elements(color) EQ 0 THEN color = 20
IF N_Elements(xsize) EQ 0 THEN xsize = 400
IF N_Elements(ysize) EQ 0 THEN ysize = 12
self.autoupdate = Keyword_Set(autoupdate)

   ; Update self structure.

self.delay = delay
self.step = theStep
self.nsteps = nsteps
self.message_label = message_label
self.message_info = message_info
self.message_xsize = message_xsize
self.message_ysize = message_ysize
self.title = title
self.color = color
self.xsize = xsize
self.ysize = ysize
self.count = 0
self.cancel = 0
self.cancelButton = Keyword_Set(cancelButton)

   ; Create the widgets.

IF N_Elements(parent) EQ 0 THEN BEGIN
   self.tlb = Widget_Base(Title=self.title, Column=1, Base_Align_Center=1)
   self.parent = -1L
ENDIF ELSE BEGIN
   IF self.cancelButton THEN modal = 0 ELSE modal = 1
   self.tlb = Widget_Base(Group_Leader=parent, Modal=modal, Title=self.title, $
      Column=1, Base_Align_Center=1, Floating=1)
   self.parent = parent
ENDELSE

self.info_baseID = Widget_Base(self.tlb, Column=1, Base_Align_Center=1)
self.acsproc_process_baseID = Widget_Base(self.tlb, Column=3)
self.message_labelID = Widget_label(self.info_baseID, Value=self.message_label, /Dynamic_Resize, /align_left,ysize=24,xsize=750,frame=1)
self.message_infoID = Widget_text(self.info_baseID, Value=self.message_info, ysize=self.message_ysize,xsize=self.message_xsize,/align_left,/scroll)
self.drawID = Widget_Draw(self.acsproc_process_baseID, XSize=self.xsize, YSize=self.ysize,/align_left)
Widget_Control, self.tlb, Set_UValue=self

IF self.cancelButton THEN BEGIN
   self.cancelID = Widget_Label(self.acsproc_process_baseID, Value='                                              ',/align_right)
   self.cancelID = Widget_Button(self.acsproc_process_baseID, Value='Cancel',/align_right)
ENDIF ELSE self.cancelID = -1L

   ; Center the top-level base.

Device, Get_Screen_Size=screenSize
xCenter = screenSize(0) / 2
yCenter = screenSize(1) / 2

geom = Widget_Info(self.tlb, /Geometry)
xHalfSize = geom.Scr_XSize / 2
yHalfSize = geom.Scr_YSize / 2

Widget_Control, self.tlb, XOffset = xCenter-xHalfSize, $
   YOffset = yCenter-yHalfSize

RETURN, 1
END; -----------------------------------------------------------------------------


PRO acsproc_process__define

; The acsproc_process class definition.

   struct = {acsproc_process, $    ; The acsproc_process object class.
             tlb:0L, $          ; The identifier of the top-level base.
             info_baseID:0L, $  ; The identifier of the info base.
             acsproc_process_baseID:0L, $  ; The identifier of the show progress base.
             button_baseID:0L, $          ; The identifier of the button base.
             message_labelID:0L, $      ; The identifier of the message_label widget.
             message_infoID:0L, $      ; The identifier of the message_label widget.
             drawID:0L, $       ; The identifier of the draw widget.
             parent:0L, $       ; The identifier of the group leader widget.
             cancelID:0L, $     ; The identifier of the CANCEL button.
             wid:0L, $          ; The window index number of the draw widget.
             xsize:0L, $        ; The XSize of the progress bar.
             ysize:0L, $        ; The YSize of the progress bar.
             color:0L, $        ; The color of the progress bar.
             autoupdate:0L, $   ; A flag for indicating if the bar should update itself.
             cancel:0L, $       ; A flag to indicate the CANCEL button was clicked.
             cancelButton:0L, $ ; A flag to indicate a CANCEL button should be added.
             message_label:'', $      ; The message_label to be written over the progress bar.
             message_info:'', $      ; The message_info to be written over the progress bar.
             message_xsize:0L, $      ; The message_info box xsize.
             message_ysize:0L, $      ; The message_info box ysize.
             title:'', $        ; The title of the top-level base widget.
             count:0L, $        ; The number of times the progress bar has been updated.
             startTime:0D, $    ; The time when the widget is started.
             percent:0D, $      ; The current percent of progress
             delay:0L, $        ; The total time the widget is on the display.
             nsteps:100L, $       ; The number of steps you want to take.
             step:0.0 $         ; The time delay between steps.
             }

END; -----------------------------------------------------------------------------
