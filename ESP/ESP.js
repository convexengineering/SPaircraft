// ESP.js implements functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer and Bob Haimes

// Copyright (C) 2010/2017  John F. Dannenhoffer, III (Syracuse University)
//
// This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//    MA  02110-1301  USA

// functions expected by wv
//    wvInitUI()                called by wv-render.js
//    wvUpdateUI()              called by wv-render.js
//    wvServerMessage(text)     called by wv-socket.js
//    wvServerDown()            called by wv-socket.js
//    wvUpdateCanvas(gl)        called by ESP.html

// functions associated with button presses (and associated button presses)
//    activateBuildButton()
//    cmdBuild()
//    cmdUndo()
//    cmdFile()
//    cmdFileNew()
//    cmdFileOpen()
//    cmdFileSave()
//    cmdFileEdit()
//       editCsmOk()
//       editCsmCancel()
//    cmdSketch()
//    cmdSketchSave()
//    cmdSketchQuit()
//    cmdStepThru()
//    cmdHelp()
//    cmdTest()
//    cmdHome()
//    cmdLeft()
//    cmdRite()
//    cmdBotm()
//    cmdTop()
//    cmdIn()
//    cmdOut()

// functions associated with menu selections (and associated button presses)
//    addPmtr()
//    editPmtr(e)
//       addRow()
//       addColumn()
//       compSens()
//       setVel(e)
//       clrVels()
//       editPmtrOk()
//       editPmtrCancel()
//    delPmtr()
//    addBrch()
//       addBrchOk()
//       addBrchCancel()
//    editBrch(e)
//       addAttr()
//       delBrch(e)
//       delAttr()
//       showBrchAttrs()
//       showBrchArgs()
//       enterSketcher()
//       buildTo()
//       editBrchOk()
//       editBrchCancel()
//    showBodyAttrs(e)

// functions associated with the mouse in the canvas
//    getMouseDown0(e)
//    getMouseMove0(e)
//    getMouseUp0(e)
//    getMouseRoll0(e)
//    mouseLeftCanvas0(e)

// function associated with the mouse in the Sketcher
//    getMouseDown7(e)
//    getMouseMove7(e)
//    getMouseUp7(e)

//  function associated with key presses in whole document
//    getKeyPress(e)
//    getKeyDown(e)
//    getKeyUp(e)

// functions associated with the mouse in the key window
//    setKeyLimits(e)

// functions associated with toggling display settings
//    toggleViz(e)
//    toggleGrd(e)
//    toggleTrn(e)
//    toggleOri(e)
//    cancelStepThru()
//    modifyDisplayType(e)
//    modifyDisplayFilter(e)
//    displayFilterOff()

// functions associated with the Sketcher
//    initializeSketch()
//    loadSketch(begs, vars, cons, segs)
//    sketchKeyPress()
//    solveSketchPre()
//    solveSketchPost(text)
//    quitSketch()
//    saveSketch()
//    undoSketch()
//    saveSketchUndo()
//
//    updateSegData(iseg)
//    drawSketch()
//    drawBezier(context, Xbezier, Ybezier)
//    drawSpline(context, Xspline, Yspline)
//    getClosestSketchPoint()
//    getClosestSketchSegment()
//    getClosestSketchConstraint()

// functions associated with a Tree in the treefrm
//    Tree(doc, treeId) - constructor
//    TreeAddNode(iparent, name, tooltip, gprim, click,
//                prop1, cbck1,
//                prop2, cbck2,
//                prop3, cbck3)
//    TreeBuild()
//    TreeClear()
//    TreeContract(inode)
//    TreeExpand(inode)
//    TreeProp(inode, iprop, onoff)
//    TreeUpdate()

// helper functions
//    resizeFrames()            called byESP.html
//    changeMode(newMode)
//    rebuildTreeWindow()
//    postMessage(mesg)
//    setupEditBrchForm()
//    setupEditPmtrForm()
//    browserToServer(text)
//    numberOfPmtrChanges()
//    numberOfBrchChanges()
//    unhighlightColumn1()
//    printObject(obj)
//    sprintf()

"use strict";


//
// callback when the user interface is to be initialized (called by wv-render.js)
//
function wvInitUI()
{
    // alert("wvInitUI()");

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix   = new J3DIMatrix4();
    wv.saveMatrix = new J3DIMatrix4(wv.mvMatrix);

                                   // ui cursor variables
    wv.cursorX   = -1;             // current cursor position
    wv.cursorY   = -1;
    wv.keyPress  = -1;             // last key pressed
    wv.keyCode   = -1;
    wv.startX    = -1;             // start of dragging position
    wv.startY    = -1;
    wv.button    = -1;             // button pressed
    wv.modifier  =  0;             // modifier (shift,alt,cntl) bitflag
    wv.flying    =  1;             // flying multiplier (do not set to 0)
    wv.offTop0   =  0;             // offset to upper-left corner of the canvas
    wv.offLeft0  =  0;
    wv.offTop7   =  0;             // offset to upper-left corner of Sketcher
    wv.offLeft7  =  0;
    wv.dragging  =  false;         // true during drag operation
    wv.picking   =  0;             // keycode of command that turned picking on
    wv.locating  =  0;             // keycode of command that turned locating on
    wv.focus     = [0, 0, 0, 1];   // focus data needed in locating
    wv.debugUI   =  0;             // set to 1 for console messages
    wv.idntStat  =  0;             // -1 server is identified
                                   //  0 need to identify server
                                   // >0 waiting for server to identify
    wv.pmtrStat  =  0;             // -2 latest Parameters are in Tree
                                   // -1 latest Parameters not in Tree (yet)
                                   //  0 need to request Parameters
                                   // >0 waiting for Parameters (request already made)
    wv.brchStat  =  0;             // -2 latest Branches are in Tree
                                   // -1 latest Branches not in Tree (yet)
                                   //  0 need to request Branches
                                   // >0 waiting for Branches (request already made)
    wv.builtTo   = 99999;          // last Branch in previous successful build
    wv.menuEvent = undefined;      // event associated with click in Tree
    wv.server    = undefined;      // string passed back from "identify;"
    wv.plotType  =  0;             // =0 mono, =1 ubar, =2 vbar, =3 cmin, =4 cmax, =5 gc
    wv.loLimit   = -1;             // lower limit in key
    wv.upLimit   = +1;             // upper limit in key
    wv.nchanges  = 0;              // number of Branch or Parameter changes by browser
    wv.filename  = "";             // name of the .csm file
    wv.curFile   = "";             // contents of current .csm file
    wv.curMode   =  0;             //-1 to disable buttons, etc.
                                   // 0 show WebViewer in canvas
                                   // 1 show addBrchForm
                                   // 2 show editBrchForm with addBrchHeader
                                   // 3 show editBrchForm with editBrchHeader
                                   // 4 show editPmtrForm with addPmtrHeader
                                   // 5 show editPmtrForm with editPmtrHeader
                                   // 6 show editCsmForm
                                   // 7 show sketcherForm
    wv.curStep   =  0;             // =1 if in StepThru mode
    wv.curPmtr   = -1;             // Parameter being editted (or -1)
    wv.curBrch   = -1;             // Branch being editted (or -1)
    wv.afterBrch = -1;             // Branch to add after (or -1)
    wv.numArgs   =  0;             // number of arguments in wv.curBrch
    wv.sgData    = {};             // scene graph metadata
    wv.scale     =  1;             // scale factor for axes
    wv.getFocus  = undefined;      // entry to get first focus
//  wv.centerV                     // set to 1 to center view and rotation
//  wv.pick                        // set to 1 to turn picking on
//  wv.locate                      // set to 1 to turn locating on
//  wv.sceneGraph                  // pointer to sceneGraph
//  wv.picked                      // sceneGraph object that was picked
//  wv.located                     // sceneGraph object that was located
//  wv.sceneUpd                    // should be set to 1 to re-render scene
//  wv.sgUpdate                    // =1 if the sceneGraph has been updated
//  wv.canvasID                    // ID of main canvas
//  wv.canvasKY                    // ID of key  canvas
//  wv.drawKey                     // =1 if key window needs to be redrawn
//  wv.socketUt.send(text)         // function to send text to server
//  wv.plotAttrs                   // plot attributes

//  sket.ibrch                     // Branch index that launched Sketcher
//  sket.mode                      // 0 initiaizing
//                                 // 1 drawing
//                                 // 2 setting curvature
//                                 // 3 constraining
//                                 // 4 setting width
//                                 // 5 setting depth
//                                 // 6 solved
//  sket.movingPoint = -1;         // moving point while constraining (or -1 for none)
//  sket.basePoint = undefined;    // base Point for W or D constraint
//  sket.skbegX    = undefined;    // expression for X from SKBEG
//  sket.skbegY    = undefined;    // expression for Y from SKBEG
//  sket.skbegZ    = undefined;    // expression for Z from SKBEG
//  sket.relative  = 0;            // =1 if relative   from SKBEG
//  sket.suggest   = undefined;    // suggested changes
//  sket.halo      =  5;           // pixels for determining "closeness" while drawing
//  sket.offTop    =  0;           // offset to upper-left corner of the canvas
//  sket.offLeft   =  0;

//  sket.scale     = undefined;    // drawing scale factor
//  sket.xorig     = undefined;    // x-origin for drawing
//  sket.yorig     = undefined;    // y-origin for drawing

//  sket.pnt.x[   i]               // x-coordinate
//  sket.pnt.y[   i]               // y-coordinate
//  sket.pnt.lbl[ i]               // constraint label

//  sket.seg.type[i]               // Segment type (L, C, S, or P)
//  sket.seg.ibeg[i]               // Point at beginning
//  sket.seg.iend[i]               // Point at end
//  sket.seg.xm[  i]               // x-coordinate at middle
//  sket.seg.ym[  i]               // y-coordinate at middle
//  sket.seg.lbl[ i]               // constraint label
//
//  sket.seg.dip[ i]               // dip (>0 for convex, <0 for concave)
//  sket.seg.xc[  i]               // x-coordinate at center
//  sket.seg.yc[  i]               // y-coordinate at center
//  sket.seg.rad[ i]               // radius (screen coordinates)
//  sket.seg.tbeg[i]               // theta at beginning (radians)
//  sket.seg.tend[i]               // theta at end       (radians)

//  sket.var.name[ i]              // variable name
//  sket.var.value[i]              // variable value

//  sket.con.type[  i]             // constraint type (X, Y, P, T, H, V, I, L, R, or S)
//  sket.con.index1[i]             // primary Point or Segment index
//  sket.con.index2[i]             // secondary Point or Segment index
//  sket.con.value[ i]             // (optional) constraint value

    sket.maxundo =  5;             // maximum   undos in Sketcher
//  sket.nundo                     // number of undos in Sketcher
//  sket.undo[i].mode              // undo snapshot (or undefined)
//  sket.undo[i].pnt
//  sket.undo[i].seg
//  sket.undo[i].var
//  sket.undo[i].con

    document.addEventListener('keypress', getKeyPress,     false);
    document.addEventListener('keydown',  getKeyDown,      false);
    document.addEventListener('keyup',    getKeyUp,        false);

    var canvas = document.getElementById(wv.canvasID);
    canvas.addEventListener('mousedown',  getMouseDown0,    false);
    canvas.addEventListener('mousemove',  getMouseMove0,    false);
    canvas.addEventListener('mouseup',    getMouseUp0,      false);
    canvas.addEventListener("wheel",      getMouseRoll0,    false);
    canvas.addEventListener('mouseout',   mouseLeftCanvas0, false);

    var sketcher = document.getElementById("sketcher");
    sketcher.addEventListener('mousedown',  getMouseDown7,  false);
    sketcher.addEventListener('mousemove',  getMouseMove7,  false);
    sketcher.addEventListener('mouseup',    getMouseUp7,    false);

    var keycan = document.getElementById(wv.canvasKY);
    keycan.addEventListener('mouseup',    setKeyLimits,    false);

    document.getElementById("sketchMenuBtn").hidden = true;
}


//
// callback when the user interface should be updated (called by wv-render.js)
//
function wvUpdateUI()
{
    // alert("wvUpdateUI()");

    // special code for delayed-picking mode
    if (wv.picking > 0) {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of 'A' operation
            if (wv.picking == 65) {
                while (1) {             // used to jump out on error

                    //  get Attribute name and value
                    var newAttrName = prompt("Enter new Attribute name:");
                    if (newAttrName === null) {
                        break;
                    }

                    var newAttrValue = prompt("Enter new Attribute value" +
                                              "(either a semi-colon separated list" +
                                              " or begin with a '$' for a string)");
                    if (newAttrValue === null) {
                        break;
                    }

                    // add statements after last executed Branch or at end of Branches
                    if (wv.builtTo < 99999) {
                        var ibrch = Number(wv.builtTo);
                    } else {
                        var ibrch = brch.length;
                    }

                    // determine the Body to select
                    var ibody = -1;
                    var gprimList = wv.picked.gprim.trim().split(" ");
                    if (gprimList.length == 3) {
                        var tempList = wv.sgData[gprimList[0]];
                        for (var ii = 0; ii < tempList.length; ii+=2) {
                            if (tempList[ii] == "_body") {
                                ibody = Number(tempList[ii+1]);
                                break;
                            }
                        }
                    } else if (gprimList.length == 4) {
                        ibody = Number(gprimList[1]);
                    }
                    if (ibody <= 0) {
                        alert("Could not find the Body");
                        break;
                    }

                    // find the _faceID or _edgeID Attributes
                    var attrs = wv.sgData[wv.picked.gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i] == "_faceID") {

                            // create the "select body" statement
                            browserToServer("newBrch|"+ibrch+"|select|$body|"+ibody+"|||||||||");

                            // create the "select face" statement
                            var faceIDlist = attrs[i+1].trim().split(" ");
                            browserToServer("newBrch|"+(ibrch+1)+"|select|$face|"+faceIDlist[0]+"|"+
                                                                                  faceIDlist[1]+"|"+
                                                                                  faceIDlist[2]+"||||||||");

                            // add the Attribute to the new select statement
                            browserToServer("setAttr|"+(ibrch+2)+"|"+newAttrName+"|1|"+newAttrValue+"|");

                            postMessage("Adding: select $body "+ibody+"\n" +
                                        "        select $face "+faceIDlist[0]+" "+
                                                                faceIDlist[1]+" "+
                                                                faceIDlist[2]+"\n" +
                                        "            attribute "+newAttrName+" "+newAttrValue+"\n" +
                                        "====> Re-build is needed <====");

                            // get an updated version of the Branches and activate Build button
                            wv.brchStat = 0;

                            activateBuildButton();
                            break;
                        } else if (attrs[i] == "_edgeID") {

                            // create the "select body" statement
                            browserToServer("newBrch|"+ibrch+"|select|$body|"+ibody+"|||||||||");

                            // create the "select edge" statement
                            var edgeIDlist = attrs[i+1].trim().split(" ");
                            browserToServer("newBrch|"+(ibrch+1)+"|select|$edge|"+edgeIDlist[0]+"|"+
                                                                                  edgeIDlist[1]+"|"+
                                                                                  edgeIDlist[2]+"|"+
                                                                                  edgeIDlist[3]+"|"+
                                                                                  edgeIDlist[4]+"||||||");

                            // add the Attribute to the new select statement
                            browserToServer("setAttr|"+(ibrch+2)+"|"+newAttrName+"|1|"+newAttrValue+"|");

                            postMessage("Adding: select $body "+ibody+"\n" +
                                        "        select $edge "+edgeIDlist[0]+" "+
                                                                edgeIDlist[1]+" "+
                                                                edgeIDlist[2]+" "+
                                                                edgeIDlist[3]+" "+
                                                                edgeIDlist[4]+"\n" +
                                        "            attribute "+newAttrName+" "+newAttrValue+"\n" +
                                        "====> Re-build is needed <====");

                            // get an updated version of the Branches and activate Build button
                            wv.brchStat = 0;

                            activateBuildButton();
                            break;
                        }
                    }
                    break;
                }

            // second part of 'g' operation
            } else if (wv.picking == 103) {
                postMessage("Toggling grid of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
                            myTree.prop(inode, 2, "on");
                        } else {
                            myTree.prop(inode, 2, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 'o' operation
            } else if (wv.picking == 111) {
                postMessage("Toggling orientation of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
                            myTree.prop(inode, 3, "on");
                        } else {
                            myTree.prop(inode, 3, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 't' operation
            } else if (wv.picking == 116) {
                postMessage("Toggling transparency of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
                            myTree.prop(inode, 3, "on");
                        } else {
                            myTree.prop(inode, 3, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 'v' operation
            } else if (wv.picking == 118) {
                postMessage("Toggling visibility of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                            myTree.prop(inode, 1, "on");
                        } else {
                            myTree.prop(inode, 1, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of '^' operation
            } else if (wv.picking == 94) {
                var mesg = "Picked: "+wv.picked.gprim;

                try {
                    var attrs = wv.sgData[wv.picked.gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
                    }
                } catch (x) {
                }
                postMessage(mesg);
            }

            wv.picked  = undefined;
            wv.picking = 0;
            wv.pick    = 0;

        // abort picking on mouse motion
        } else if (wv.dragging) {
            postMessage("Picking aborted");

            wv.picking = 0;
            wv.pick    = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }

    // special code for delayed-locating mode
    if (wv.locating > 0) {

        // if something is located, post a message
        if (wv.located !== undefined) {

            // second part of '@' operation
            if (wv.locating == 64) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Located: x="+xloc.toFixed(4)+", y="+yloc.toFixed(4)+
                                   ", z="+zloc.toFixed(4));
            }

            wv.located  = undefined;
            wv.locating = 0;
            wv.locate   = 0;

        // abort locating on mouse motion
        } else if (wv.dragging) {
            postMessage("Locating aborted");

            wv.locating = 0;
            wv.locate   = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }

    // if the server is not identified, try to identify it now
    if (wv.idntStat > 0) {
        wv.idntStat--;
    } else if (wv.idntStat == 0) {
        try {
            browserToServer("identify|");
            browserToServer("getFilename|");
            wv.idntStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.idntStat = 10;
       }
    }

    // if the Parameters are scheduled to be updated, send a message to
    //    get the Parameters now
    if (wv.pmtrStat > 0) {
        wv.pmtrStat--;
    } else if (wv.pmtrStat == 0) {
        try {
            browserToServer("getPmtrs|");
            wv.pmtrStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.pmtrStat = 10;
        }
    }

    // if the Branches are scheduled to be updated, send a message to
    //    get the Branches now
    if (wv.brchStat > 0) {
        wv.brchStat--;
    } else if (wv.brchStat == 0) {
        try {
            browserToServer("getBrchs|");
            wv.brchStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.brchStat = 10;
        }
    }

    // if the scene graph and Parameters have been updated, (re-)build the Tree
    if ((wv.sgUpdate == 1 && wv.pmtrStat <= -1 && wv.brchStat <= -1) ||
        (                    wv.pmtrStat == -1 && wv.brchStat == -2) ||
        (                    wv.pmtrStat == -2 && wv.brchStat == -1)   ) {

        if (wv.sceneGraph === undefined) {
            alert("wv.sceneGraph is undefined --- but we need it");
        }

        rebuildTreeWindow();
    }

    // deal with key presses
    if (wv.keyPress != -1 && wv.curMode == 0) {

        var myKeyPress = String.fromCharCode(wv.keyPress);

        // '?' -- help
        if (myKeyPress == "?") {
            postMessage("........................... Viewer Cursor options ...........................\n" +
                        "ctrl-h <Home> - initial view             ctrl-f        - front view          \n" +
                        "ctrl-l        - leftside view            ctrl-r        - riteside view       \n" +
                        "ctrl-t        - top view                 ctrl-b        - bottom view         \n" +
                        "ctrl-i <PgUp> - zoom in                  ctrl-o <PgDn> - zoom out            \n" +
                        "<Left>        - rotate or xlate left     <Rite>        - rotate or xlate rite\n" +
                        "<Up>          - rotate or xlate up       <Down>        - rotate or xlate down\n" +
                        ">             - save view                <             - recall view         \n" +
                        "ctrl->        - save view to file        ctrl-<        - read view from file \n" +
                        "^             - query object at cursor   @             - get coords. @ cursor\n" +
                        "v             - toggle Viz at cursor     g             - toggle Grd at cursor\n" +
                        "t             - toggle Trn at cursor     o             - toggle Ori at cursor\n" +
                        "A             - add Attribute at cursor  *             - center view @ cursor\n" +
                        "!             - toggle flying mode       ?             - get help            \n" +
                        ".............................................................................");

        // 'A' -- add Attribute at cursor
        } else if (myKeyPress == 'A') {
            wv.picking  = 65;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'g' -- toggle grid at cursor
        } else if (myKeyPress == 'g') {
            wv.picking  = 103;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'o' -- orientation at cursor
        } else if (myKeyPress == 'o') {
            wv.picking  = 111;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 't' -- toggle transparency at cursor
        } else if (myKeyPress == 't') {
            wv.picking  = 116;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'v' -- toggle visibility at cursor
        } else if (myKeyPress == 'v') {
            wv.picking  = 118;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '^' -- query at cursor
        } else if (myKeyPress == '^') {
            wv.picking  = 94;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '@' -- locate at cursor
        } else if (myKeyPress == '@') {
            wv.locating = 64;
            wv.locate   = 1;

        // '*' -- center view
        } else if (myKeyPress == '*') {
            wv.centerV = 1;
            postMessage("View being centered");

        // '!' -- toggle flying mode
        } else if (myKeyPress == "!") {
            if (wv.flying <= 1) {
                postMessage("Turning flying mode ON");
                wv.flying = 10;
            } else {
                postMessage("Turning flying mode OFF");
                wv.flying = 1;
            }

        // '>' -- save view
        } else if (myKeyPress == ">") {
            postMessage("Saving current view");
            wv.saveMatrix.load(wv.mvMatrix);
            wv.sceneUpd = 1;

        // '<' -- recall view
        } else if (myKeyPress == "<") {
            postMessage("Restoring saved view");
            wv.mvMatrix.load(wv.saveMatrix);
            wv.sceneUpd = 1;

        // C-'>' -- save view to file
        } else if (wv.keyPress == 46 && wv.modifier == 5) {
            var filename = prompt("Enter view filename:");
            if (filename !== null) {
                postMessage("Saving view to \"" + filename + "\"");
                browserToServer("saveView|"+filename+"|"+wv.scale+"|"+wv.mvMatrix.getAsArray()+"|");
            }

        // C-'<' -- read view from file
        } else if (wv.keyPress == 44 && wv.modifier == 5) {
            var filename = prompt("enter view filename:");
            if (filename !== null) {
                postMessage("Reading view from \"" + filename + "\"");
                browserToServer("readView|"+filename+"|");
            }

        // '<esc>' -- not used
//      } else if (wv.keyPress == 0 && wv.keyCode == 27) {
//          postMessage("<Esc> is not supported.  Use '?' for help");

        // '<Home>' -- initial view
        } else if (wv.keyPress == 0 && wv.keyCode == 36) {
            wv.mvMatrix.makeIdentity();
            wv.scale    = 1;
            wv.sceneUpd = 1;

        // '<end>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 35) {
            postMessage("<End> is not supported.  Use '?' for help");

        // '<PgUp>' -- zoom in
        } else if (wv.keyPress == 0 && wv.keyCode == 33) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(2.0, 2.0, 2.0);
                wv.scale *= 2.0;
            } else {
                wv.mvMatrix.scale(1.25, 1.25, 1.25);
                wv.scale *= 1.25;
            }
            wv.sceneUpd = 1;

        // '<PgDn>' -- zoom out
        } else if (wv.keyPress == 0 && wv.keyCode == 34) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(0.5, 0.5, 0.5);
                wv.scale *= 0.5;
            } else {
                wv.mvMatrix.scale(0.8, 0.8, 0.8);
                wv.scale *= 0.8;
            }
            wv.sceneUpd = 1;

        // '<Delete>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 46) {
            postMessage("<Delete> is not supported.  Use '?' for help");

        // '<Left>' -- rotate or translate object left
        } else if (wv.keyPress == 0 && wv.keyCode == 37) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( -5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(-0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(-0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Right>' -- rotate or translate object right
        } else if (wv.keyPress == 0 && wv.keyCode == 39) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( +5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(+0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(+0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Up>' -- rotate or translate object up
        } else if (wv.keyPress == 0 && wv.keyCode == 38) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( -5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, +0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, +0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Down>' -- rotate or translate object down
        } else if (wv.keyPress == 0 && wv.keyCode == 40) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( +5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, -0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, -0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // 'ctrl-h' - initial view (same as <Home>)
        } else if ((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode ==  8)   ) {
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-i' - zoom in (same as <PgUp> without shift)
        } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
            cmdIn();
            wv.keyPress = -1;
            return;

        // '+' - zoom in (same as <PgUp> without shift)
        } else if (wv.keyPress ==  43 && wv.modifier == 1) {
            cmdIn();
            wv.keyPress = -1;
            return;

        // 'ctrl-o' - zoom out (same as <PgDn> without shift)
        } else if ((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) {
            cmdOut();
            wv.keyPress = -1;
            return;

        // '-' - zoom out (same as <PgDn> without shift)
        } else if (wv.keyPress ==  45 && wv.modifier == 0) {
            cmdOut();
            wv.keyPress = -1;
            return;

        // 'ctrl-f' - front view (same as <Home>)
        } else if ((wv.keyPress == 102 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   6 && wv.modifier == 4 && wv.keyCode ==  6)   ) {
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-r' - riteside view
        } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
            cmdRite();
            wv.keyPress = -1;
            return;

        // 'ctrl-l' - leftside view
        } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
            cmdLeft();
            wv.keyPress = -1;
            return;

        // 'ctrl-t' - top view
        } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
            cmdTop();
            wv.keyPress = -1;
            return;

        // 'ctrl-b' - bottom view
        } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
            cmdBotm();
            wv.keyPress = -1;
            return;

            // NOP
        } else if (wv.keyPress == 0 && wv.modifier == 0) {

        // unknown command
//      } else {
//          postMessage("Unknown command (keyPress="+wv.keyPress
//                      +", modifier="+wv.modifier
//                      +", keyCode="+wv.keyCode+").  Use '?' for help");
        }

        wv.keyPress = -1;
    }

    // UI is in screen coordinates (not object)
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();

    // deal with mouse movement
    if (wv.dragging) {

        // cntrl is down (rotate)
        if (wv.modifier == 4) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
            if ((angleX != 0.0) || (angleY != 0.0)) {
                wv.mvMatrix.rotate(angleX, 1,0,0);
                wv.mvMatrix.rotate(angleY, 0,1,0);
                wv.sceneUpd = 1;
            }

        // alt-shift is down (rotate)
        } else if (wv.modifier == 3) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
            if ((angleX != 0.0) || (angleY != 0.0)) {
                wv.mvMatrix.rotate(angleX, 1,0,0);
                wv.mvMatrix.rotate(angleY, 0,1,0);
                wv.sceneUpd = 1;
            }

        // alt is down (spin)
        } else if (wv.modifier == 2) {
            var xf = wv.startX - wv.width  / 2;
            var yf = wv.startY - wv.height / 2;

            if ((xf != 0.0) || (yf != 0.0)) {
                var theta1 = Math.atan2(yf, xf);
                xf = wv.cursorX - wv.width  / 2;
                yf = wv.cursorY - wv.height / 2;

                if ((xf != 0.0) || (yf != 0.0)) {
                    var dtheta = Math.atan2(yf, xf) - theta1;
                    if (Math.abs(dtheta) < 1.5708) {
                        var angleZ = 128*(dtheta) / 3.1415926 / wv.flying;
                        wv.mvMatrix.rotate(angleZ, 0,0,1);
                        wv.sceneUpd = 1;
                    }
                }
            }

        // shift is down (zoom)
        } else if (wv.modifier == 1) {
            if (wv.cursorY != wv.startY) {
                var scale = Math.exp((wv.cursorY - wv.startY) / 512.0 / wv.flying);
                wv.mvMatrix.scale(scale, scale, scale);
                wv.scale   *= scale;
                wv.sceneUpd = 1;
            }

        // no modifier (translate)
        } else {
            var transX = (wv.cursorX - wv.startX) / 256.0 / wv.flying;
            var transY = (wv.cursorY - wv.startY) / 256.0 / wv.flying;
            if ((transX != 0.0) || (transY != 0.0)) {
                wv.mvMatrix.translate(transX, transY, 0.0);
                wv.sceneUpd = 1;
            }
        }

        // if not flying, then update the start coordinates
        if (wv.flying <= 1) {
            wv.startX = wv.cursorX;
            wv.startY = wv.cursorY;
        }
    }
}


//
// callback when a (text) message is received from the server (called by wv-socket.js)
//
function wvServerMessage(text)
{
    // alert("wvServerMessage(text="+text+")");

    if (wv.debugUI) {
        var date = new Date;
        console.log("("+date.toTimeString().substring(0,8)+") browser<--server: "+text.substring(0,40));
    }

    // if it is a null message, do nothing
    if (text.length <= 1) {

    // if it starts with "identify|" post a message */
    } else if (text.substring(0,9) == "identify|") {
        if (wv.server === undefined) {
            wv.server = text.substring(9,text.length-2);
            postMessage("ESP has been initialized and is attached to '"+wv.server+"'");
        }

    // if it starts with "sgData|" store the auxiliary scene graph data
    } else if (text.substring(0,7) == "sgData|") {
        if (text.length > 9) {

            // convert  ".]" to "]"  and  ",}" to "}"
            var newtext = text.substring(7,text.length-1);
            newtext = newtext.replace(/,]/g, "]").replace(/,}/g, "}");
            wv.sgData = JSON.parse(newtext);
        } else {
            wv.sgData = {};
        }

    // if it starts with "sgFocus|" store the scene graph focus data
    } else if (text.substring(0,8) == "sgFocus|") {
        if (text.length > 10) {
            wv.focus = JSON.parse(text.substring(8,text.length-1));
        } else {
            wv.focus = [];
        }

    // if it starts with "getPmtrs|" build the (global) pmtr array
    } else if (text.substring(0,9) == "getPmtrs|") {
        if (text.length > 11) {
            pmtr = JSON.parse(text.substring(9,text.length-1));
        } else {
            pmtr = new Array;
        }
        wv.pmtrStat = -1;

        rebuildTreeWindow();

    // if it starts with "newPmtr|" do nothing
    } else if (text.substring(0,8) == "newPmtr|") {

    // if it starts with "clrVels|" do nothing
    } else if (text.substring(0,8) == "clrVels|") {

    // if it starts with "setPmtr|" do nothing
    } else if (text.substring(0,8) == "setPmtr|") {

    // if it starts with "delPmtr|" do nothing
    } else if (text.substring(0,8) == "delPmtr|") {

    // if it starts with "setVel|" do nothing
    } else if (text.substring(0,7) == "setVel|") {

    // if it starts with "getBrchs|" build the (global) brch array
    } else if (text.substring(0,9) == "getBrchs|") {
        if (text.length > 11) {
            brch = JSON.parse(text.substring(9,text.length-1));
        } else {
            brch = new Array;
        }
        wv.brchStat = -1;

        rebuildTreeWindow();

    // if it starts with "newBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "newBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8,text.length-1));
        }

    // if it starts with "setBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "setBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8,text.length-1));
        }

    // if it starts with "delBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "delBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8,text.length-1));
        }

    // if it starts with "setAttr|" do nothing
    } else if (text.substring(0,8) == "setAttr|") {

    // if it starts with "undo|" get the latest Parameters and
    //                           Branches and update tree (scene)
    } else if (text.substring(0,5) == "undo|") {
        var cmd = text.substring(5,text.length-2);

        wv.sgUpdate = 1;
        wv.brchStat = 0;
        wv.pmtrStat = 0;

        postMessage("Undoing '"+cmd+"' ====> Re-build is needed <====");

        activateBuildButton();

    // if it starts with "new|" post message
    } else if (text.substring(0,4) == "new|") {
        postMessage("New (blank) configuration has been loaded\n" +
                    "    Press \"Design Parameters\" (in left window) to add a Design Parameter\n" +
                    "    Press \"Branches\" (in left window) to begin a 3D Object\n" +
                    "    Press \"Branches\" and then \"skbeg\" to begin a 2D Sketch");

        pmtr = new Array();
        brch = new Array();

        rebuildTreeWindow();

    // if it starts with "save|" do nothing
    } else if (text.substring(0,5) == "save|") {

    // if it starts with "build|" reset the build button
    } else if (text.substring(0,6) == "build|") {
        var button = document.getElementById("buildButton");
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;

        var textList = text.split("|");

        if (textList[1].substring(0,7) == "ERROR::") {
            postMessage(textList[1]);

            var ibrch    = Number(textList[2]);
            var nbody    = Number(textList[3]);

            wv.builtTo = ibrch;

            if (nbody > 0) {
                postMessage("Build complete through ibrch="+ibrch+
                            ", which generated "+nbody+" Body(s)");
            } else {
                alert("No Bodys were produced");
            }
        } else {
            var ibrch    = Number(textList[1]);
            var nbody    = Number(textList[2]);

            wv.builtTo = ibrch;

            if (ibrch == brch.length) {
                postMessage("Entire build complete, which generated "+nbody+
                            " Body(s)");
            } else if (ibrch >= brch.length) {
                postMessage("Build complete through ibrch="+ibrch+
                            ", which generated "+nbody+" Body(s)");
            } else if (ibrch > 0) {
                postMessage("Partial build (through "+brch[ibrch-1].name+
                            ") complete, which generated "+nbody+" Body(s)");
            } else {
                postMessage("textList[1]="+textList[1]);
                postMessage("ibrch="+ibrch+"   brch.length="+brch.length);
                postMessage("Build failed in first Branch");
            }
        }

        changeMode(0);

    // if it starts with "loadSketch|" initialize the Sketcher
    } else if (text.substring(0,11) == "loadSketch|") {

        var textList = text.split("|");

        if (textList[1].substring(0,7) == "ERROR::") {
            postMessage(textList[1]);
            alert("Sketch cannot be loaded");
        } else {
            // open the Sketcher
            changeMode(7);
            initializeSketch();

            loadSketch(textList[1], textList[2], textList[3], textList[4]);
        }

    // if it starts with "solveSketch|" update the Sketcher
    } else if (text.substring(0,12) == "solveSketch|") {

        solveSketchPost(text);

    // if it starts with "saveSketchBeg|" do nothing
    } else if (text.substring(0,14) == "saveSketchBeg|") {

    // if it starts with "saveSketchMid|" do nothing
    } else if (text.substring(0,14) == "saveSketchMid|") {

    // if it starts with "saveSketchEnd|" do nothing
    } else if (text.substring(0,14) == "saveSketchEnd|") {

    // if it starts with "saveSketch|" either return an error or update the mode
    } else if (text.substring(0,11) == "saveSketch|") {

        if (text.substring(11,14) == "ok|") {
            activateBuildButton();
            changeMode(0);
        } else {
            alert("The sketch could not be saved");
        }

    // if it starts with "setLims|" do nothing
    } else if (text.substring(0,8) == "setLims|") {

    // if it starts with "getFilename|" store the results in wv.filename
    } else if (text.substring(0,12) == "getFilename|") {
        var textList = text.split("|");

        wv.filename = textList[1];

        if (wv.filename.length <= 0) {
            postMessage("ESP has started without a .csm file\n" +
                        "    Press \"Design Parameters\" (in left window) to add a Design Parameter\n" +
                        "    Press \"Branches\" (in left window) to begin a 3D Object\n" +
                        "    Press \"Branches\" and then \"skbeg\" to begin a 2D Sketch");
        } else {
            postMessage("\"" + wv.filename + "\" has been loaded");
        }

    // if it starts with "getCsmFile|" store the results in wv.curFile and change to mode 6
    } else if (text.substring(0,11) == "getCsmFile|") {
        var ichar = 11;

        wv.curFile = text.substring(ichar);

        if (wv.filename.length <= 0) {
            postMessage("ESP has started without a .csm file\n" +
                        "    Press \"Design Parameters\" (in left window) to add a Design Parameter\n" +
                        "    Press \"Branches\" (in left window) to begin a 3D Object\n" +
                        "    Press \"Branches\" and then \"skbeg\" to begin a 2D Sketch");
        }

        // fill in the name of the .csm file
        var csmFilename = document.getElementById("editCsmFilename");

        csmFilename["innerHTML"] = "Contents of: "+wv.filename;

        // fill the textarea with the current .csm file
        var csmTextArea = document.getElementById("editCsmTextArea");

        csmTextArea.cols  = 84;
        csmTextArea.rows  = 25;
        csmTextArea.value = wv.curFile;

        // post the editCsmForm
        changeMode(6);

    // if it starts with "setCsmFileBeg|" do nothing
    } else if (text.substring(0,14) == "setCsmFileBeg|") {

    // if it starts with "setCsmFileMid|" do nothing
    } else if (text.substring(0,14) == "setCsmFileMid|") {

    // if it starts with "setCsmFileEnd|" do nothing
    } else if (text.substring(0,14) == "setCsmFileEnd|") {

    // if it starts with "setWvKey|" turn key or logo on
    } else if (text.substring(0,9) == "setWvKey|") {
        if (wv.curMode == 0) {
            if (text.substring(9,11) == "on") {
                document.getElementById("WVkey"  ).hidden = false;
                document.getElementById("ESPlogo").hidden = true;
            } else {
                document.getElementById("WVkey"  ).hidden = true;
                document.getElementById("ESPlogo").hidden = false;
            }
        }

    // if it starts with "saveView|" do nothing
    } else if (text.substring(0,9) == "saveView|") {

    // if it starts with "readView|" load view matrix and update display
    } else if (text.substring(0,9) == "readView|") {
        var readViewList = text.split("|");
        var entries      = readViewList[2].split(",");
        if (entries.length == 16) {
            var matrix  = Array(16);
            for (var i = 0; i < 16; i++) {
                matrix[i] = parseFloat(entries[i]);
            }
            wv.mvMatrix.makeIdentity();
            wv.uiMatrix.load(matrix);
            wv.scale    = parseFloat(readViewList[1]);
            wv.sceneUpd = 1;
        } else {
            postMessage("File has wrong number of entries");
        }

    // if it starts with "nextStep|" either switch button legend back
    //  to StepThru or post message about what is showing
    } else if (text.substring(0,9) == "nextStep|") {
        var nextStepList = text.split("|");
        if (nextStepList[2] == "") {
            wv.curStep = 0;
            var button = document.getElementById("stepThruBtn");
            button["innerHTML"] = "StepThru";

            postMessage("Finished with StepThru mode");
        } else {
            wv.curStep = 1;
            postMessage("Showing \""+nextStepList[1]+"\" generated by \""+
                       nextStepList[2]+"\" in StepThru mode");
        }

    // if it starts with "ERROR:: there is nothing to undo" post the error
    } else if (text.substring(0,32) == "ERROR:: there is nothing to undo") {
        alert("There is nothing to undo");

    // if it starts with "ERROR:: save(" post the error
    } else if (text.substring(0,13) == "ERROR:: save(") {
        postMessage(text);
        alert("File not written");

    // if it starts with "ERROR::" post the error
    } else if (text.substring(0,7) == "ERROR::") {
        postMessage(text);
        alert("Error encountered (see messages)");

        var button = document.getElementById("buildButton");

        button["innerHTML"] = "Fix before re-build";
        button.style.backgroundColor = "#FF3F3F";

        changeMode(0);

    // default is to post the message
    } else {
        postMessage(text);
    }
}


//
// callback when the server goes down either on abort or normal closure (called by wv-socket.js)
//
function wvServerDown()
{
    // deactivate the buttons
    wv.curMode = -1;

    postMessage("The server has terminated or network connection has been lost.\n"
                +"ESP will be working off-line");

    // turn the background of the message window pink
    var botm = document.getElementById("brframe");
    botm.style.backgroundColor = "#FFAFAF";
}


//
// callback used to put axes on the canvas (called by ESP.html)
//
function wvUpdateCanvas(gl)
{
    // Construct the identity as projection matrix and pass it in
    wv.mvpMatrix.load(wv.perspectiveMatrix);
    wv.mvpMatrix.setUniform(gl, wv.u_modelViewProjMatrixLoc, false);

    var mv   = new J3DIMatrix4();
    var mVal = wv.mvMatrix.getAsArray();

    mVal[ 3] = 0.0;
    mVal[ 7] = 0.0;
    mVal[11] = 0.0;
    mv.load(mVal);
    mv.scale(1.0/wv.scale, 1.0/wv.scale, 1.0/wv.scale);
    mv.invert();
    mv.transpose();

    // define location of axes in space
    var x    = -1.5 * wv.width / wv.height;
    var y    = -1.5;
    var z    =  0.9;
    var alen =  0.25;     // length of axes

    // set up coordinates for axes
    mVal = mv.getAsArray();

    var vertices = new Float32Array(66);
    vertices[ 0] = x;
    vertices[ 1] = y;
    vertices[ 2] = z;
    vertices[ 3] = x + alen*(    mVal[ 0]             );
    vertices[ 4] = y + alen*(    mVal[ 1]             );
    vertices[ 5] = z + alen*(    mVal[ 2]             );
    vertices[ 6] = x + alen*(1.1*mVal[ 0]+0.1*mVal[ 4]);
    vertices[ 7] = y + alen*(1.1*mVal[ 1]+0.1*mVal[ 5]);
    vertices[ 8] = z + alen*(1.1*mVal[ 2]+0.1*mVal[ 6]);
    vertices[ 9] = x + alen*(1.3*mVal[ 0]-0.1*mVal[ 4]);
    vertices[10] = y + alen*(1.3*mVal[ 1]-0.1*mVal[ 5]);
    vertices[11] = z + alen*(1.3*mVal[ 2]-0.1*mVal[ 6]);
    vertices[12] = x + alen*(1.1*mVal[ 0]-0.1*mVal[ 4]);
    vertices[13] = y + alen*(1.1*mVal[ 1]-0.1*mVal[ 5]);
    vertices[14] = z + alen*(1.1*mVal[ 2]-0.1*mVal[ 6]);
    vertices[15] = x + alen*(1.3*mVal[ 0]+0.1*mVal[ 4]);
    vertices[16] = y + alen*(1.3*mVal[ 1]+0.1*mVal[ 5]);
    vertices[17] = z + alen*(1.3*mVal[ 2]+0.1*mVal[ 6]);

    vertices[18] = x;
    vertices[19] = y;
    vertices[20] = z;
    vertices[21] = x + alen*(    mVal[ 4]             );
    vertices[22] = y + alen*(    mVal[ 5]             );
    vertices[23] = z + alen*(    mVal[ 6]             );
    vertices[24] = x + alen*(1.1*mVal[ 4]+0.1*mVal[ 8]);
    vertices[25] = y + alen*(1.1*mVal[ 5]+0.1*mVal[ 9]);
    vertices[26] = z + alen*(1.1*mVal[ 6]+0.1*mVal[10]);
    vertices[27] = x + alen*(1.2*mVal[ 4]             );
    vertices[28] = y + alen*(1.2*mVal[ 5]             );
    vertices[29] = z + alen*(1.2*mVal[ 6]             );
    vertices[30] = x + alen*(1.3*mVal[ 4]+0.1*mVal[ 8]);
    vertices[31] = y + alen*(1.3*mVal[ 5]+0.1*mVal[ 9]);
    vertices[32] = z + alen*(1.3*mVal[ 6]+0.1*mVal[10]);
    vertices[33] = x + alen*(1.2*mVal[ 4]             );
    vertices[34] = y + alen*(1.2*mVal[ 5]             );
    vertices[35] = z + alen*(1.2*mVal[ 6]             );
    vertices[36] = x + alen*(1.2*mVal[ 4]             );
    vertices[37] = y + alen*(1.2*mVal[ 5]             );
    vertices[38] = z + alen*(1.2*mVal[ 6]             );
    vertices[39] = x + alen*(1.2*mVal[ 4]-0.1*mVal[ 8]);
    vertices[40] = y + alen*(1.2*mVal[ 5]-0.1*mVal[ 9]);
    vertices[41] = z + alen*(1.2*mVal[ 6]-0.1*mVal[10]);


    vertices[42] = x;
    vertices[43] = y;
    vertices[44] = z;
    vertices[45] = x + alen*(    mVal[ 8]             );
    vertices[46] = y + alen*(    mVal[ 9]             );
    vertices[47] = z + alen*(    mVal[10]             );
    vertices[48] = x + alen*(1.1*mVal[ 8]+0.1*mVal[ 0]);
    vertices[49] = y + alen*(1.1*mVal[ 9]+0.1*mVal[ 1]);
    vertices[50] = z + alen*(1.1*mVal[10]+0.1*mVal[ 2]);
    vertices[51] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
    vertices[52] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
    vertices[53] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
    vertices[54] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
    vertices[55] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
    vertices[56] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
    vertices[57] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
    vertices[58] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
    vertices[59] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
    vertices[60] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
    vertices[61] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
    vertices[62] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
    vertices[63] = x + alen*(1.3*mVal[ 8]-0.1*mVal[ 0]);
    vertices[64] = y + alen*(1.3*mVal[ 9]-0.1*mVal[ 1]);
    vertices[65] = z + alen*(1.3*mVal[10]-0.1*mVal[ 2]);

    // set up colors for the axes
    var colors = new Uint8Array(66);
    colors[ 0] = 255;   colors[ 1] =   0;   colors[ 2] =   0;
    colors[ 3] = 255;   colors[ 4] =   0;   colors[ 5] =   0;
    colors[ 6] = 255;   colors[ 7] =   0;   colors[ 8] =   0;
    colors[ 9] = 255;   colors[10] =   0;   colors[11] =   0;
    colors[12] = 255;   colors[13] =   0;   colors[14] =   0;
    colors[15] = 255;   colors[16] =   0;   colors[17] =   0;
    colors[18] =   0;   colors[19] = 255;   colors[20] =   0;
    colors[21] =   0;   colors[22] = 255;   colors[23] =   0;
    colors[24] =   0;   colors[25] = 255;   colors[26] =   0;
    colors[27] =   0;   colors[28] = 255;   colors[29] =   0;
    colors[30] =   0;   colors[31] = 255;   colors[32] =   0;
    colors[33] =   0;   colors[34] = 255;   colors[35] =   0;
    colors[36] =   0;   colors[37] = 255;   colors[38] =   0;
    colors[39] =   0;   colors[40] = 255;   colors[41] =   0;
    colors[42] =   0;   colors[43] =   0;   colors[44] = 255;
    colors[45] =   0;   colors[46] =   0;   colors[47] = 255;
    colors[48] =   0;   colors[49] =   0;   colors[50] = 255;
    colors[51] =   0;   colors[52] =   0;   colors[53] = 255;
    colors[54] =   0;   colors[55] =   0;   colors[56] = 255;
    colors[57] =   0;   colors[58] =   0;   colors[59] = 255;
    colors[60] =   0;   colors[61] =   0;   colors[62] = 255;
    colors[63] =   0;   colors[64] =   0;   colors[65] = 255;

    // draw the axes
    if (gl.lineWidth) {
        gl.lineWidth(2);
    }
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.LINES, 0, 22);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);
}


//
// activate the buildButton
//
function activateBuildButton() {
    // alert("in activateBuildButton()");

    var button = document.getElementById("buildButton");

    button["innerHTML"] = "Press to Re-build";
    button.style.backgroundColor = "#3FFF3F";
}


//
// callback when "buildButton" is pressed
//
function cmdBuild() {
    // alert("in cmdBuild()");

    var button  = document.getElementById("buildButton");
    var buttext = button["innerHTML"];

    if (buttext == "Up to date") {
        if (confirm("The configuration is up to date.\n" +
                    "Do you want to force a rebuild?") === true) {
            postMessage("Forced re-building...");

            // build first so that parameters are updated
            browserToServer("build|-1|");

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;

            browserToServer("getBrchs|");
            wv.brchStat = 6000;

            button["innerHTML"] = "Re-building...";
            button.style.backgroundColor = "#FFFF3F";

            //inactivate buttons until build is done
            changeMode(-1);
        }

    } else if (buttext == "Press to Re-build") {
        if (wv.curMode != 0) {
            alert("Command disabled,  Press 'Cancel' or 'OK' first");
            return;
        }

        // build first so that parameters are updated
        browserToServer("build|0|");

        browserToServer("getPmtrs|");
        wv.pmtrStat = 6000;

        browserToServer("getBrchs|");
        wv.brchStat = 6000;

        button["innerHTML"] = "Re-building...";
        button.style.backgroundColor = "#FFFF3F";

        //inactivate buttons until build is done
        changeMode(-1);

    } else if (buttext == "Re-building...") {
        alert("Rebuilding in process.  Please be patient");

    } else if (buttext == "Fix before re-build") {
        alert("Edit/add/delete a Branch to fix before re-building");

    } else if (buttext == "Press to Solve") {
        if (wv.curMode != 7) {
            alert("Command disabled,  Press 'Cancel' or 'OK' first");
            return;
        }

        // save undo info
        saveSketchUndo();

        // send solve request to the server
        solveSketchPre();

    } else if (buttext == "Initializing...") {
        alert("The Sketcher is initializing.  Please be patient");

    } else if (buttext == "Drawing...") {
        alert("Use \"L\" to add a line segment\n" +
              "    \"C\" to add a circular arc\n" +
              "    \"S\" to add a spline point\n" +
              "    \"B\" to add a Bezier point\n" +
              "    \"Z\" to add a zero-length segment, or\n" +
              "    \"O\" to leave the Sketch open");

    } else if (buttext == "Setting R...") {
        alert("Press left mouse to set radius");

    } else if (buttext == "Constraining...") {
        // save undo info
        saveSketchUndo();

        // send solve request to the server
        solveSketchPre();

    } else if (buttext == "nothing") {
        alert("If at a point, use\n" +
              "    \"X\" to set the X-coordinate\n" +
              "    \"Y\" to set the Y-coordinate\n" +
              "    \"P\" to set perpendicularity\n" +
              "    \"T\" to set tangency\n" +
              "    \"A\" to set angle (positive to left)\n" +
              "    \"W\" to set width (dx) to another point\n" +
              "    \"D\" to set depth (dy) to another point\n" +
              "    \"<\" to delete a constraint\n" +
              "If at a cirarc segment, use\n" +
              "    \"R\" to set radius\n" +
              "    \"S\" to set sweep ccw sweep angle\n" +
              "    \"X\" to set X-coordinate of center\n" +
              "    \"Y\" to set Y-coordinate of center\n" +
              "    \"<\" to delete\n" +
              "If at a line segment, use\n" +
              "    \"H\" to set a horizontal constraint\n" +
              "    \"V\" to set a vertical constraint\n", +
              "    \"I\" to set inclination (from right horizontal)\n" +
              "    \"L\" to set length\n"+
              "    \"<\" to delete\n" +
              "If at a W or D constraint, use\n" +
              "    \"<\" to delete");

    } else if (buttext == "Setting W...") {
        alert("Move to another point and press left mouse to set width (dx)");

    } else if (buttext == "Setting D...") {
        alert("move to another point and press left mouse to set depth (dy)");

    } else {
        alert("Unexpected button text:\""+buttext+"\"");
    }
}


//
// callback when "undoButton" is pressed (called by ESP.html)
//
function cmdUndo() {
    // alert("in cmdUndo()");

    if (wv.curMode == 7) {
        undoSketch();
        return;
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("cmdUndo is not implemented for "+wv.server);
        return;
    }

    browserToServer("undo|");
}


//
// callback when "fileButton" is pressed (called by ESP.html)
//
function cmdFile() {
    // alert("in cmdFile()");

    // toggle between hiding and showing the File menu contents
    document.getElementById("myFileMenu").classList.toggle("showFileMenu");
}


//
// callback when "File->New" is pressed (called by ESP.html)
//
function cmdFileNew() {
    // alert("in cmdFileNew()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.server != "serveCSM") {
        alert("cmdFileNew is not implemented for "+wv.server);
        return;

    } else if (wv.nchanges > 0) {
        if (confirm(wv.nchanges+" change(s) will be lost.  Continue?") !== true) {
            return;
        }
    }

    if (wv.curMode == 0) {
        browserToServer("new|");

        wv.filename = "";
        wv.nchanges = 0;

        pmtr   = new Array();
        brch   = new Array();
        sgData = {};
    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
}


//
// callback when "File->Open" is pressed (called by ESP.html)
//
function cmdFileOpen() {
    // alert("in cmdFileOpen()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.server != "serveCSM") {
        alert("cmdFileOpen is not implemented for "+wv.server);
        return;

    } else if (wv.nchanges > 0) {
        if (confirm(wv.nchanges+" change(s) will be lost.  Continue?") !== true) {
            return;
        }
    }

    if (wv.curMode == 0) {
        var filename = prompt("Enter filename:", wv.filename);
        if (filename !== null) {
            if (filename.search(/\.csm$/) > 0 ||
                filename.search(/\.cpc$/) > 0 ||
                filename.search(/\.udc$/) > 0   ) {
                // well-formed filename
            } else {
                // add .csm extension
                filename += ".csm";
            }

            postMessage("Attempting to open \""+filename+"\" ...");

            browserToServer("open|"+filename+"|");

            wv.filename = filename;
            wv.nchanges = 0;

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;

            browserToServer("getBrchs|");
            wv.brchStat = 6000;

            var button = document.getElementById("buildButton");
            button["innerHTML"] = "Re-building...";
            button.style.backgroundColor = "#FFFF3F";

            //inactivate buttons until build is done
            changeMode(-1);
        } else {
            postMessage("NOT opening since no filename specified");
        }

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
}


//
// callback when "File->Save" is pressed (called by ESP.html)
//
function cmdFileSave() {
    // alert("in cmdFileSave()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.server != "serveCSM") {
        alert("cmdFileSave is not implemented for "+wv.server);
        return;

    } else if (wv.curMode == 0) {
        var filename = prompt("Enter filename:", wv.filename);
        if (filename !== null) {
            if (filename.search(/\.csm$/) > 0 ||
                filename.search(/\.cpc$/) > 0 ||
                filename.search(/\.udc$/) > 0   ) {
                // well-formed filename
            } else {
                // add .csm extension
                filename += ".csm";
            }

            postMessage("Saving model to '"+filename+"'");
            browserToServer("save|"+filename+"|");

            wv.filename = filename;
            wv.nchanges = 0;
        } else {
            postMessage("NOT saving since no filename specified");
        }

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
}


//
// callback when "editButton" is pressed (called by ESP.html)
//
function cmdFileEdit() {
    // alert("in cmdFileEdit()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("cmdFileEdit is not implemented for "+wv.server);
        return;
    }

    // if there have been any changes, tell user to save first
    if (wv.nchanges > 0) {
        alert("Changes have been made to the 'Design Parameters'\n" +
              "and/or the 'Branches' via the user-interface.\n\n" +
              "Save these changes before you can edit the .csm file");
        return;
    }

    // get the latest csm file
    try {
        browserToServer("getCsmFile|");
        wv.nchanges = 0;
    } catch (e) {
        // could not send message
    }
}


//
// callback when "OK" button is pressed in editCsmForm (called by ESP.html)
//
function editCsmOk() {
    // alert("in editCsmOk()");

    var editCsmTextArea = document.getElementById("editCsmTextArea");

    // tell user if no changes were made
    if (wv.curFile == editCsmTextArea.value) {
        alert("No changes were made");
        changeMode(0);
        return;
    } else {
        if (wv.filename.length == 0) {
            wv.filename = prompt("Enter filename");
            if (wv.filename !== null) {
                if (wv.filename.search(/\.csm$/) <= 0) {
                    // add .csm extension
                    wv.filename += ".csm";
                }
            } else {
                alert("NOT saving since no filename specified");
                return;
            }
        } else if (confirm("This will overwrite your input file.  Continue?") !== true) {
            return;
        }
    }

    // because of an apparent limit on the size of text
    //    messages that can be sent from the browser to the
    //    server, we need to send the new file back in
    //    pieces and then reassemble on the server
    var maxMessageSize = 800;

    var ichar = 0;
    var part  = editCsmTextArea.value.substring(ichar, ichar+maxMessageSize);
    browserToServer("setCsmFileBeg|"+wv.filename+"|"+part);
    ichar += maxMessageSize;

    while (ichar < editCsmTextArea.value.length) {
        part = editCsmTextArea.value.substring(ichar, ichar+maxMessageSize);
        browserToServer("setCsmFileMid|"+part);
        ichar += maxMessageSize;
    }

    browserToServer("setCsmFileEnd|");

    // remember the edited .csm file
    wv.curFile = editCsmTextArea.value;

    postMessage("'"+wv.filename+"' file has been changed.");

    // get an updated version of the Parameters and Branches
    wv.pmtrStat = 0;
    wv.brchStat = 0;

    // remove the contents of the file from memory
    wv.curFile = "";

    // reset the number of changes
    wv.nchanges = 0;

    // inform the user that a rebuild is in process
    var button = document.getElementById("buildButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode( 0);
    changeMode(-1);
}


//
// callback when "Cancel" is pressed in editCsmForm (called by ESP.html)
//
function editCsmCancel() {
    // alert("in editCsmCancel()");

    // remove the contents of the file from memory
    wv.curFile = "";

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "sketchButton" is pressed (called by ESP.html)
//
function cmdSketch() {
    // alert("in cmdSketch()");

    // toggle between hiding and showing the File menu contents
    document.getElementById("mySketchMenu").classList.toggle("showSketchMenu");
}


//
// callback when "Sketch->Save" is pressed (called by ESP.html)
//
function cmdSketchSave() {
    // alert("in cmdSketchSave()");

    // close the Sketch menu
    var menu = document.getElementsByClassName("sketchMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showSketchMenu")) {
            menu[i].classList.remove(  "showSketchMenu");
        }
    }

    if (wv.server != "serveCSM") {
        alert("cmdSave is not implemented for "+wv.server);
        return;

    } else if (wv.curMode == 7) {
        saveSketch();
        wv.brchStat = 0;

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
}


//
// callback when "Sketch->Quit" is pressed (called by ESP.html)
//
function cmdSketchQuit() {
    // alert("in cmdSketchQuit()");

    // close the Sketch menu
    var menu = document.getElementsByClassName("sketchMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showSketchMenu")) {
            menu[i].classList.remove(  "showSketchMenu");
        }
    }

    if (wv.curMode != 7) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("cmdQuit is not implemented for "+wv.server);
        return;
    }

    quitSketch();
    activateBuildButton();
}


//
// callback when "stepThruBtn" is pressed (called by ESP.html)
//
function cmdStepThru() {
    // alert("in cmdStepThru()");

    if (wv.curMode >= 0) {
        browserToServer("nextStep|");

        var button = document.getElementById("stepThruBtn");
        button["innerHTML"] = "NextStep";
    } else {
        alert("Button disabled");
    }
}


//
// callback when "helpButton" is pressed (called by ESP.html)
//
function cmdHelp() {

    // open help in another tab
    window.open("ESP-help.html");
}


//
// callback when "testButton" is pressed (called by ESP.html)
//
function cmdTest() {
    alert("in cmdTest()")

    browserToServer("animate|");
}


//
// callback when "homeButton" is pressed (called by ESP.html)
//
function cmdHome()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        // get the extrema of the data
        var xmin = +1.0e+10;
        var xmax = -1.0e+10;
        var ymin = +1.0e+10;
        var ymax = -1.0e+10;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            xmin = Math.min(xmin, sket.scale * (sket.pnt.x[ipnt] - sket.xorig));
            xmax = Math.max(xmax, sket.scale * (sket.pnt.x[ipnt] - sket.xorig));
            ymin = Math.min(ymin, sket.scale * (sket.yorig - sket.pnt.y[ipnt]));
            ymax = Math.max(ymax, sket.scale * (sket.yorig - sket.pnt.y[ipnt]));
        }

        // get the size of the canvas
        var canvas = document.getElementById("sketcher");
        var width  = canvas.clientWidth;
        var height = canvas.clientHeight;

        var test1 = (xmax - xmin) / width;
        var test2 = (ymax - ymin) / height;

        // set up sizing so that Sketch fills middle 50 percent of the window
        var newScale = 2 * Math.max(test1, test2);
        var newXorig = (width  - (xmin + xmax) / newScale) / 2;
        var newYorig = (height + (ymin + ymax) / newScale) / 2;

        // convert the Points by the new scale factors
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            var xx = sket.scale * (sket.pnt.x[ipnt] - sket.xorig);
            var yy = sket.scale * (sket.yorig - sket.pnt.y[ipnt]);

            sket.pnt.x[ipnt] = Math.floor(newXorig + xx / newScale);
            sket.pnt.y[ipnt] = Math.floor(newYorig - yy / newScale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            sket.seg.dip[iseg] = sket.seg.dip[iseg] * sket.scale / newScale;

            updateSegData(iseg);
        }

        sket.scale = newScale;
        sket.xorig = newXorig;
        sket.yorig = newYorig;

        drawSketch();
    }
}


//
// callback when "leftButton" is pressed (called by ESP.html)
//
function cmdLeft()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(+90, 0,1,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newXorig = sket.xorig - 100;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var xx = sket.scale * (sket.pnt.x[ipnt] - sket.xorig);
            sket.pnt.x[ipnt] = Math.floor(newXorig + xx / sket.scale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            updateSegData(iseg);
        }

        sket.xorig = newXorig;

        drawSketch();
    }
}


//
// callback when "riteButton" is pressed (called by ESP.html)
//
function cmdRite()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(-90, 0,1,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newXorig = sket.xorig + 100;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var xx = sket.scale * (sket.pnt.x[ipnt] - sket.xorig);
            sket.pnt.x[ipnt] = Math.floor(newXorig + xx / sket.scale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            updateSegData(iseg);
        }

        sket.xorig = newXorig;

        drawSketch();
    }
}


//
// callback when "botmButton" is pressed (called by ESP.html)
//
function cmdBotm()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(-90, 1,0,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newYorig = sket.yorig + 100;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var yy = sket.scale * (sket.yorig - sket.pnt.y[ipnt]);
            sket.pnt.y[ipnt] = Math.floor(newYorig - yy / sket.scale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            updateSegData(iseg);
        }

        sket.yorig = newYorig;

        drawSketch();
    }
}


//
// callback when "topButton" is pressed (called by ESP.html)
//
function cmdTop()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(+90, 1,0,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newYorig = sket.yorig - 100;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var yy = sket.scale * (sket.yorig - sket.pnt.y[ipnt]);
            sket.pnt.y[ipnt] = Math.floor(newYorig - yy / sket.scale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            updateSegData(iseg);
        }

        sket.yorig = newYorig;

        drawSketch();
    }
}


//
// callback when "inButton" is pressed (called by ESP.html)
//
function cmdIn()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.load(wv.uiMatrix);
        wv.mvMatrix.scale(2.0, 2.0, 2.0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale   *= 2.0;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newScale = 0.5 * sket.scale;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var xx = sket.scale * (sket.pnt.x[ipnt] - sket.xorig);
            var yy = sket.scale * (sket.yorig - sket.pnt.y[ipnt]);

            sket.pnt.x[ipnt] = Math.floor(sket.xorig + xx / newScale);
            sket.pnt.y[ipnt] = Math.floor(sket.yorig - yy / newScale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            sket.seg.dip[iseg] = sket.seg.dip[iseg] * sket.scale / newScale;

            updateSegData(iseg);
        }

        sket.scale = newScale;

        drawSketch();
    }
}


//
// callback when "outButton" is pressed (called by ESP.html)
//
function cmdOut()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.load(wv.uiMatrix);
        wv.mvMatrix.scale(0.5, 0.5, 0.5);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale   *= 0.5;
        wv.sceneUpd = 1;
    } else if (wv.curMode == 7) {
        var npnt = sket.pnt.x.length;
        var nseg = sket.seg.type.length;

        if (npnt <= 1) {
            alert("Not enough data to rescale");
            return;
        } else if (sket.scale === undefined || isNaN(sket.scale)) {
            alert("At least one length or radius must be set before rescaling  (sket.scale="+sket.scale+")");
            return;
        } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
            alert("At lest one \"X\" must be set before rescaling  (sket.xorig="+sket.xorig+")");
            return;
        } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
            alert("At lest one \"Y\" must be set before rescaling  (sket.yorig="+sket.yorig+")");
            return;
        }

        var newScale = 2.0 * sket.scale;

        for (var ipnt = 0; ipnt < npnt; ipnt++) {
            var xx = sket.scale * (sket.pnt.x[ipnt] - sket.xorig);
            var yy = sket.scale * (sket.yorig - sket.pnt.y[ipnt]);

            sket.pnt.x[ipnt] = Math.floor(sket.xorig + xx / newScale);
            sket.pnt.y[ipnt] = Math.floor(sket.yorig - yy / newScale);
        }

        for (var iseg = 0; iseg < nseg; iseg++) {
            sket.seg.dip[iseg] = sket.seg.dip[iseg] * sket.scale / newScale;

            updateSegData(iseg);
        }

        sket.scale = newScale;

        drawSketch();
    }
}


//
// callback when "Design Parameters" is pressed in Tree
//
function addPmtr() {
    // alert("addPmtr()");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("addPmtr is not implemented for "+wv.server);
        return;
    }

    // get the new name
    var name = prompt("Enter new Parameter name");
    if (name === null) {
        return;
    } else if (name.length <= 0) {
        return;
    }

    // check that name is valid
    if (name.match(/^[a-zA-Z]\w*$/) === null) {
        alert("'"+name+"' is not a valid name");
        return;
    }

    // check that the name does not exist already
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (name == pmtr[ipmtr].name) {
            alert("'"+name+"' already exists");
            return;
        }
    }

    // store the values locally
    var newPmtr = pmtr.length;

    pmtr[newPmtr] = new Array();

    pmtr[newPmtr].name = name;
    pmtr[newPmtr].type = 500;
    pmtr[newPmtr].nrow = 1;
    pmtr[newPmtr].ncol = 1;
    pmtr[newPmtr].value = new Array(1);

    pmtr[newPmtr].value[0] = "";

    // remember info for Parameter
    wv.curPmtr = newPmtr;

    // set up editPmtrForm
    if (setupEditPmtrForm() > 0) {
        return;
    }

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
}


//
// callback when Design Parameter name is pressed in Tree
//
function editPmtr(e) {
    // alert("in editPmtr(e="+e+")");

    if        (wv.curMode == 5) {
        // currently editting another Parameter, so cancel (throwing away changes)
        editPmtrCancel();
    } else if (wv.curMode == 3) {
        // currently editting a Branch,          so cancel (throwing away changes)
        editBrchCancel();
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("editPmtr is not implemented for "+wv.server);
        return;
    }

    wv.menuEvent = e;

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }
    if (ipmtr < 0) {
        alert(name+" not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // highlight this Parameter in the Tree
    var myElem = document.getElementById("node"+inode+"col1");
    myElem.className = "currentTD";

    // highlight all Branches that explicitly mention this Parameter
    var re = RegExp("[\\][)(+\\-*/^.,;]"+pmtr[ipmtr].name+"[\\][)(+\\-*/^.,;]");
    for (var ibrch = 0; ibrch < brch.length; ibrch++) {
        for (var iarg = 0; iarg < brch[ibrch].args.length; iarg++) {
            if (("("+brch[ibrch].args[iarg]+")").match(re) !== null) {
                for (var jnode = 0; jnode < myTree.name.length; jnode++) {
                    if (myTree.parent[jnode] == 3 &&
                        myTree.name[jnode].replace(/\u00a0/g, "") == brch[ibrch].name) {
                        document.getElementById("node"+jnode+"col1").className = "childTD";
                    }
                }
            }
        }
    }

    // remember info for the current Parameter
    wv.curPmtr = ipmtr;

    // set up editPmtrForm
    setupEditPmtrForm();

    // post the editPmtr form (with the editPmtr header)
    changeMode(5);
}


//
// callback when "Add row" is pressed in editPmtrForm (called by ESP.html)
//
function addRow() {
    // alert("in addRow()");

    // adjust the number of rows
    pmtr[wv.curPmtr].nrow++;
    pmtr[wv.curPmtr].value = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);

    for (var i = 0; i < pmtr[wv.curPmtr].value.length; i++) {
        pmtr[wv.curPmtr].value[i] = "";
    }

    // set up editPmtrForm
    if (setupEditPmtrForm() > 0) {
        return;
    }

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
}


//
// callback when "Add column" is pressed in editPmtrForm (called by ESP.html)
//
function addColumn() {
    // alert("in addColumn()");

    // adjust the number of columns
    pmtr[wv.curPmtr].ncol++;
    pmtr[wv.curPmtr].value = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);

    for (var i = 0; i < pmtr[wv.curPmtr].value.length; i++) {
        pmtr[wv.curPmtr].value[i] = "";
    }

    // set up editPmtrForm
    if (setupEditPmtrForm() > 0) {
        return;
    }

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
}


//
// callback when "Compute Sensitivity" is pressed in editPmtrForm (called by ESP.html)
//
function compSens() {
    // alert("in compSens()");

    // disable this command if there were any changes to the Parameter
    if (numberOfPmtrChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }
    if (ipmtr < 0) {
        alert(name+" not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // can only compute sensitivity for a scalar
    if (pmtr[ipmtr].nrow > 1 || pmtr[ipmtr].ncol > 1) {
        alert("Use \"Set Design Velocity\" to select which element of this multi-valued parameter to use.  Then \"Press to Re-build\"");
        return;
    }

    // clear any previous velocities
    browserToServer("clrVels|");
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        pmtr[jpmtr].dot[0] = 0;
    }

    // set velocity for ipmtr
    pmtr[ipmtr].dot[0] = 1.0;
    browserToServer("setVel|"+pmtr[ipmtr].name+"|1|1|1|");
    postMessage("Computing sensitivity with respect to "+pmtr[ipmtr].name);

    // rebuild
    browserToServer("build|0|");
    wv.brchStat = 6000;

    browserToServer("getPmtrs|");
    wv.pmtrStat = 6000;

    browserToServer("getBrchs|");

    var button = document.getElementById("buildButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode( 0);
    changeMode(-1);
}


//
// callback when "Set Design Velocity" is pressed in editPmtrForm (called by ESP.html)
//
function setVel() {
    // alert("in setVel()");

    // disable this command if there were any changes to the Parameter
    if (numberOfPmtrChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }
    if (ipmtr < 0) {
        alert(name+" not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // get each of the values
    var index  = -1;
    var nchange = 0;
    for (var irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
        for (var icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
            index++;

            // get the new value
            var newVel;
            if (pmtr[ipmtr].nrow == 1 && pmtr[ipmtr].ncol == 1) {
                newVel = prompt("Enter new Design Velocity for "+name,
                                pmtr[ipmtr].dot[index]);
            } else {
                newVel = prompt("Enter new Design Velocity for "+name+
                                "["+irow+","+icol+"]",
                                pmtr[ipmtr].dot[index]);
            }

            // make sure a valid number was entered
            if (newVel === null) {
                continue;
            } else if (isNaN(newVel)) {
                alert("Illegal number format, so Design Velocity not being changed");
                continue;
            } else if (newVel == pmtr[ipmtr].dot[index]) {
                continue;
            } else if (pmtr[ipmtr].nrow == 1 && pmtr[ipmtr].ncol == 1) {
                postMessage("Parameter '"+name+"' has new Design Velocity "+
                            newVel+" ====> Re-build is needed <====");
                nchange++;
            } else {
                postMessage("Parameter '"+name+"["+irow+","+icol+
                            "]' has new Design Velocity "+newVel+
                            " ====> Re-build is needed <====");
                nchange++;
            }

            // store the value locally
            pmtr[ipmtr].dot[index] = Number(newVel);

            // send the new Design Velocity to the server
            browserToServer("setVel|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+newVel+"|");

        }
    }

    // update the UI
    if (nchange > 0) {
        wv.nchanges += nchange;

        var myElem = document.getElementById(id);
        myElem.className = "fakelinkoff";

        activateBuildButton();
    }

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "Clear Design Velocities" is pressed in editPmtrForm (called by ESP.html)
//
function clrVels() {
    // alert("clrVels()");

    // disable this command if there were any changes to the Parameter
    if (numberOfPmtrChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("clrVels is not implemented for "+wv.server);
        return;
    }

    // get an updated Parameter list (so that added Parameter is listed)
    browserToServer("clrVels|");

    // update the UI
    postMessage("Design Velocities have been cleared ====> Re-build is needed <====");
    activateBuildButton();
}


//
// callback when "OK" is pressed in editPmtrForm (called by ESP.html)
//
function editPmtrOk() {
    // alert("in editPmtrOk()");

    var editPmtrForm = document.getElementById("editPmtrForm");

    var ipmtr = wv.curPmtr;
    var name  = pmtr[ipmtr].name;
    var nrow  = pmtr[ipmtr].nrow;
    var ncol  = pmtr[ipmtr].ncol;
    var irow;
    var icol;

    // make sure that all entries have valid values
    var nchange = 0;
    var index   = -1;
    for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
        for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
            index++;

            // get the new value
            var myInput = editPmtrForm["row"+irow+"col"+icol];
            var value   = myInput.value.replace(/\s/g, "");

            if (value.length <= 0) {
                alert("Entry at (row "+irow+", col "+icol+") is blank");
                return;
            } else if (isNaN(value)) {
                alert("Illegal number format at (row "+irow+", col "+icol+")");
                return;
            }
        }
    }

    // send the new Parameter to the server if in add Brch mode
    if (wv.curMode == 4) {
        var mesg = "newPmtr|"+name+"|"+nrow+"|"+ncol+"|";

        index = -1;
        for (irow = 1; irow <= nrow; irow++) {
            for (icol = 1; icol <= ncol; icol++) {
                index++;
                mesg = mesg+"|";
            }
        }

        browserToServer(mesg);
    }

    // get each of the values
    index = -1;
    for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
        for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
            index++;

            // get the new value
            var myInput = editPmtrForm["row"+irow+"col"+icol];
            var value = myInput.value.replace(/\s/g, "");

            if (value != pmtr[ipmtr].value[index]) {
                postMessage("Parameter '"+pmtr[ipmtr].name+"["+irow+","+icol+
                            "]' has been changed to "+value+
                            " ====> Re-build is needed <====");
                nchange++;

                // store the value locally
                pmtr[ipmtr].value[index] = Number(value);

                // send the new value to the server
                browserToServer("setPmtr|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+value+"|");
            }
        }
    }

    // update the UI
    if (nchange > 0) {
        wv.nchanges += nchange;

        if (wv.curMode != 4) {
            var id     = wv.menuEvent["target"].id;
            var myElem = document.getElementById(id);
            myElem.className = "fakelinkoff";

        // get an updated Parameter list (so that added Pmtr is listed)
        } else {
            browserToServer("getPmtrs|");
        }

        activateBuildButton();
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "Cancel" is pressed in editPmtrForm (called by ESP.html)
//
function editPmtrCancel() {
    // alert("in editPmtrCancel()");

    // if we are in process of adding a Parameter, remove it now
    if (wv.curMode == 4) {
        pmtr.splice(pmtr.length-1, 1);
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "Delete Parameter" is pressed in editPmtrForm (called by ESP.html)
//
function delPmtr() {
    // alert("in delPmtr()");

    var ipmtr = wv.curPmtr + 1;

    // send message to the server
    browserToServer("delPmtr|"+pmtr[wv.curPmtr].name+"|");

    // get updated Parameters
    browserToServer("getPmtrs|");
    wv.pmtrStat = 0;

    // update the UI
    postMessage("Deleting Parameter "+name+" ====> Re-build is needed <====");
    activateBuildButton();

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "Branch" is pressed in Tree
// callback when "Add new Branch after this Brch" is pressed in editBrchForm (called by ESP.html)
//
function addBrch() {
    // alert("in addBrch()");

    // this check allows one to select AddBranchAfterThisBranch
    if (wv.menuEvent !== undefined && wv.curMode == 3) {
        if (brch[wv.curBrch].type == "udprim") {
            var arg0 = brch[wv.curBrch].args[0];
            if (arg0.charAt(1) == "$" || arg0.charAt(1) == "/") {
                alert("Cannot add a Branch within a UDC");
                return;
            }
        }
    } else if (wv.curMode != 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("addBrch is not implemented for "+wv.server);
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // remember
    if (wv.curBrch >= 0) {
        wv.afterBrch = wv.curBrch;
        wv.curBrch   = -1;
    } else {
        wv.afterBrch = brch.length - 1;
        wv.curBrch   = -1;
    }

    // post the addBrch form
    changeMode(1);
}


//
// callback when "OK" is pressed in addBrchForm (called by ESP.html)
//
function addBrchOk() {
    // alert("in addBrchOk()");

    // get the elements on the form
    var elements = document.getElementById("addBrchForm").elements;

    var jelem       = -1;
    for (var ielem = 0; ielem < elements.length; ielem++) {
        if (elements[ielem].checked) {
            jelem = ielem;
            break;
        }
    }

    // if nothing is selected, return to the WebViewer
    if (jelem < 0) {
        alert("Select a Branch type or press 'Cancel'");
        return;
    }

    // if something was picked, initialize the new Branch
    var newBrch = brch.length;

    brch[newBrch] = new Array();

    brch[newBrch].name = "**name_automatically_assigned**";
    brch[newBrch].type = elements[jelem].value;
    brch[newBrch].actv = 300;
    brch[newBrch].args = new Array();

    // remember info for Branch to add after
    wv.curBrch = newBrch;

    // set up editBrchForm (and abort if problem is detected)
    if (setupEditBrchForm() > 0) {
        return;
    }

    // increment the number of changes
    wv.nchanges++;

    // post the editBrch form (with the addBrch header)
    changeMode(2);
}


//
// callback when "Cancel" is pressed in addBrchForm (called by ESP.html)
//
function addBrchCancel() {
    // alert("in addBrchCancel()");

    // return to the WebViewer
    changeMode(0);
}


//
// callback when Branch name is pressed in Tree
//
function editBrch(e) {
    // alert("in editBrch(e="+e+")");

    if        (wv.curMode == 3) {
        // currently editting another Branch, so cancel (throwing away changes)
        editBrchCancel();
    } else if (wv.curMode == 5) {
        // currently editting a Parameter,    so cancel (throwing away changes)
        editPmtrCancel();
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("editBrch is not implemented for "+wv.server);
        return;
    }

    changeMode(0);

    wv.menuEvent = e;

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    var name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    var ibrch = -1;           // 0-bias
    var jbrch;                // 1-bias (and temp)
    for (jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch;
            break;
        }
    }
    if (ibrch < 0) {
        alert(name+" not found");
        return;
    } else {
        jbrch = ibrch + 1;
    }

    // highlight this Branch in the Tree
    var myElem = document.getElementById("node"+inode+"col1");
    myElem.className = "currentTD";

    // highlight the parents and child of this Branch in the Tree */
    var ileft = brch[ibrch].ileft;
    var irite = brch[ibrch].irite;
    var ichld = brch[ibrch].ichld;

    for (var jnode = 0; jnode < myTree.name.length; jnode++) {
        if (ileft > 0) {
            if (myTree.parent[jnode] == 3 &&
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ileft-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "parentTD";
            }
        }
        if (irite > 0) {
            if (myTree.parent[jnode] == 3 &&
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[irite-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "parentTD";
            }
        }
        if (ichld > 0) {
            if (myTree.parent[jnode] == 3 &&
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ichld-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "childTD";
            }
        }
    }

    // remember info for the current Branch
    wv.curBrch   = ibrch;

    // set up editBrchForm
    if (setupEditBrchForm() > 0) {
        alert("This Branch cannot be editted\n(associated with a UDC)");
        unhighlightColumn1();
        changeMode(0);
        return;
    }

    // post the editBrch form (with the editBrch header)
    changeMode(3);
}


//
// callback when "Add Attribute" is pressed in editBrchForm  (called by ESP.html)
//
function addAttr() {
    // alert("in addAttr()");

    // disable this command if there were any chnges to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    var name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    var ibrch = -1;           // 0-bias
    var jbrch;                // 1-bias (and temp)
    for (jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch;
            break;
        }
    }
    if (ibrch < 0) {
        alert(name+" not found");
        return;
    } else {
        jbrch = ibrch + 1;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // get the Attribute name and value
    var atype  = prompt("Enter 1 for Attribute or 2 for Csystem", "1");
    if (atype === null) {
        return;
    } else if (atype == "1") {
        var aname  = prompt("Enter Attribute name");
        if (aname === null) {
            return;
        } else if (aname.length <= 0) {
            return;
        }
    } else if (atype == "2") {
        var aname  = prompt("Enter Csystem name");
        if (aname === null) {
            return;
        } else if (aname.length <= 0) {
            return;
        }
    } else {
        return;
    }
    var avalue = prompt("Enter Attribute value");
    if (avalue === null) {
        return;
    } else if (avalue.length <= 0) {
        return;
    }

    // send the new value to the server
    var mesg = "setAttr|"+jbrch+"|"+aname+"|"+atype+"|"+avalue+"|";

    browserToServer(mesg);

    // update the UI
    postMessage("Adding attribute '"+aname+"' (with value "+avalue+") to "+
                name+" ====> Re-build is needed <====");
    activateBuildButton();

    // get an updated version of the Branches
    wv.brchStat = 0;

    // return to  the WebViewer
    changeMode(0);
}


//
// callback when "Delete this Branch" is pressed in editBrchForm (called by ESP.html)
//
function delBrch() {
    // alert("in delBrch()");

    // disable this command if there were any changes to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    var ibrch = -1;
    var name  = "";

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if this was called by pressing "Delete last Branch", adjust inode
    //    to point to last visible node (as if "Delete this Branch"
    //    was pressed on the last Branch)
    if (inode == 2) {
        inode = myTree.child[inode];

        while (myTree.next[inode] >= 0) {
            inode = myTree.next[inode];
        }
    }

    // get the Branch name
    name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    for (var jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch + 1;
        }
    }
    if (ibrch <= 0) {
        alert(name+" not found");
        return;
    }

    // do not allow a UDPARG or UDPRIM that points to a UDC to be deleted
    if (brch[ibrch-1].type == "udparg" || brch[ibrch-1].type == "udprim") {
        if (brch[ibrch-1].args[0].charAt(1) == "/" ||
            brch[ibrch-1].args[0].charAt(1) == "$"   ) {
            alert("Cannot delete \""+brch[ibrch-1].name+"\" that points to a UDC");
            return;
        }
    }

    // send the message to the server
    browserToServer("delBrch|"+ibrch+"|");

    // hide the Tree node that was just updated
    var element = myTree.document.getElementById("node"+inode);
    element.style.display = "none";

    // get an updated version of the Branches
    wv.brchStat = 0;

    // update the UI
    postMessage("Deleting Branch "+name+" ====> Re-build is needed <====");
    activateBuildButton();

    // return to the WebViewer
    changeMode(0);
}


//
// callback when "Delete an Attribute" is pressed in editBrchForm (called by ESP.html)
//
function delAttr() {
    //alert("delAttr");

    alert("To delete an Attribute, make its value blank");
}


//
// callback when "Show Attributes" is pressed in editBrchForm (called by ESP.html)
//
function showBrchAttrs() {
    // alert("showBrchAttrs()");

    document.getElementById("AddBrchOrAttr").value   = "Add Attribute/Csystem";
    document.getElementById("AddBrchOrAttr").onclick = addAttr;

    document.getElementById("DelBrchOrAttr").value   = "Delete an Attribute/Csystem";
    document.getElementById("DelBrchOrAttr").onclick = delAttr;

    document.getElementById("ShowArgOrAttr").value   = "Show Arguments";
    document.getElementById("ShowArgOrAttr").onclick = showBrchArgs;

    document.getElementById("editBrchHeader2").innerHTML = "<h3>... or edit the attributes/csystems of the current Branch</h3>";

    document.getElementById("editBrchArgs" ).hidden = true;
    document.getElementById("editBrchAttrs").hidden = false;
}


//
// callback when "Show Arguments" is pressed in editBrchForm (called by ESP.html)
//
function showBrchArgs() {
    // alert("showBrchArgs()");

    document.getElementById("AddBrchOrAttr").value   = "Add new Branch after this Branch";
    document.getElementById("AddBrchOrAttr").onclick = addBrch;

    document.getElementById("DelBrchOrAttr").value   = "Delete this Branch";
    document.getElementById("DelBrchOrAttr").onclick = delBrch;

    document.getElementById("ShowArgOrAttr").value   = "Show Attributes/Csystems";
    document.getElementById("ShowArgOrAttr").onclick = showBrchAttrs;

    document.getElementById("editBrchHeader2").innerHTML = "<h3>... or edit the arguments of the current Branch</h3>";

    document.getElementById("editBrchArgs" ).hidden = false;
    document.getElementById("editBrchAttrs").hidden = true;
}


//
// callback when "Enter Sketcher" is pressed in editBrchForm (called by ESP.html)
//
function enterSketcher() {
    // alert("enterSketcher()");

    sket.ibrch = -1;
    var name   = "";

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    for (var jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            sket.ibrch = jbrch + 1;
        }
    }
    if (sket.ibrch <= 0) {
        alert(name+" not found");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // send the message to the server
    browserToServer("loadSketch|"+sket.ibrch+"|");

    // inactivate buttons until build is done
    changeMode(-1);

    // initialize the Sketcher
    initializeSketch();
}


//
// callback when "Build to the Branch" is pressed in editBrchForm (called by ESP.html)
//
function buildTo() {
    // alert("buildTo()");

    // disable this command if there were any chnges to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    var ibrch = -1;
    var name  = "";

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    for (var jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch + 1;
        }
    }
    if (ibrch <= 0) {
        alert(name+" not found");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // send the message to the server
    browserToServer("build|"+ibrch+"|");

    browserToServer("getPmtrs|");
    wv.pmtrStat = 6000;

    browserToServer("getBrchs|");
    wv.brchStat = 6000;

    // update the UI
    postMessage("Re-building only to "+name+"...");

    var button = document.getElementById("buildButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode(-1);
}


//
// callback when "OK" button is pressed in editBrchForm (called by ESP.html)
//
function editBrchOk() {
    // alert("in editBrchOk()");

    var editBrchForm = document.getElementById("editBrchForm");

    var ibrch = wv.curBrch;
    var jbrch = ibrch + 1;

    var brchName = editBrchForm.brchName.value.replace(/\s/g, "");
    if (brchName.length <= 0) {
        alert("Name cannot be blank");
        return;
    }

    var mesg;
    var output;
    var nchange = 0;

    if (wv.curMode == 2) {
        mesg = "newBrch|"+(wv.afterBrch+1)+"|"+brch[ibrch].type+"|";
    } else {
        mesg = "setBrch|"+jbrch+"|";
    }

    // make sure that name does not contain a space or ">"
    var newBrchName = editBrchForm.brchName.value;

    if (newBrchName.indexOf(" ") >= 0) {
        alert("Changed name '"+newBrchName+"' cannot contain a space");
        return;
    } else if (newBrchName.indexOf(">") >= 0) {
        alert("Changed name '"+newBrchName+"' cannot contain '>'");
        return;
    }

    // make sure that we are not adding or changing a Branch associated with a UDC
    if (brch[ibrch].type == "udparg" || brch[ibrch].type == "udprim") {
        var value = editBrchForm.argValu1.value.replace(/\s/g, "");
        if (value.charAt(0) == "/" || value.charAt(0) == "$") {
            if (wv.curMode == 2) {
                alert("Cannot add a \""+brch[ibrch].type+"\" that calls a UDC");
                return;
            } else if (brch[ibrch].args[0] != "$"+value) {
                alert("Cannot change primtype to a UDC");
                return;
            }
        }
    }

    if (newBrchName != brch[ibrch].name) {
        // make sure that name does not start with "Brch_"
        if (newBrchName.substring(0,5) == "Brch_") {
            alert("Changed name '"+newBrchName+
                  "' cannot begin with 'Brch_'");
            return;
        }

        // make sure that name does not already exist
        for (var kbrch = 0; kbrch < brch.length; kbrch++) {
            if (brch[kbrch].name == newBrchName) {
                alert("Name '"+newBrchName+"' already exists");
                return;
            }
        }

        // store the value locally
        postMessage("Changing name '"+brch[ibrch].name+"' to '"+
                    newBrchName+"' ====> Re-build is needed <====");
        brch[ibrch].name = newBrchName;

        nchange++;
    }

    if (wv.curMode == 3) {
        mesg = mesg + newBrchName + "|";
    }

    // update the activity
    var newActivity = editBrchForm.activity.value;
    if (brch[ibrch].actv != "none"     &&
        brch[ibrch].actv != newActivity  ) {

        // store the value locally
        brch[ibrch].actv = newActivity;
        if (brch[ibrch].actv == 300) {
            postMessage("Activating Branch "+newBrchName+
                        " ====> Re-build is needed <====");
        } else {
            postMessage("Suppressing Branch "+newBrchName+
                        " ====> Re-build is needed <====");
        }

        // get an updated version of the Branches
        wv.brchStat = 0;
        nchange++;
    }

    if (wv.curMode == 2) {

    } else if (newActivity == 300) {
        mesg = mesg + "active|";
    } else {
        mesg = mesg + "suppressed|";
    }

    // update any necessary arguments
    var prev_value = "dum";
    for (var iarg = 0; iarg < wv.numArgs; iarg++) {
        if        (iarg == 0) {
            var name  = document.getElementById("argName1").firstChild["data"];
            var value = editBrchForm.            argValu1.value.replace(/\s/g, "");
        } else if (iarg == 1) {
            var name  = document.getElementById("argName2").firstChild["data"];
            var value = editBrchForm.            argValu2.value.replace(/\s/g, "");
        } else if (iarg == 2) {
            var name  = document.getElementById("argName3").firstChild["data"];
            var value = editBrchForm.            argValu3.value.replace(/\s/g, "");
        } else if (iarg == 3) {
            var name  = document.getElementById("argName4").firstChild["data"];
            var value = editBrchForm.            argValu4.value.replace(/\s/g, "");
        } else if (iarg == 4) {
            var name  = document.getElementById("argName5").firstChild["data"];
            var value = editBrchForm.            argValu5.value.replace(/\s/g, "");
        } else if (iarg == 5) {
            var name  = document.getElementById("argName6").firstChild["data"];
            var value = editBrchForm.            argValu6.value.replace(/\s/g, "");
        } else if (iarg == 6) {
            var name  = document.getElementById("argName7").firstChild["data"];
            var value = editBrchForm.            argValu7.value.replace(/\s/g, "");
        } else if (iarg == 7) {
            var name  = document.getElementById("argName8").firstChild["data"];
            var value = editBrchForm.            argValu8.value.replace(/\s/g, "");
        } else if (iarg == 8) {
            var name  = document.getElementById("argName9").firstChild["data"];
            var value = editBrchForm.            argValu9.value.replace(/\s/g, "");
        }

        // make sure non-blank does not follow blank
        if (value.length > 0 && prev_value.length == 0) {
            alert(name+" should be blank (follows blank)");
            return;
        }

        // check for blanks (allowed for some Branch types because they
        // have variable number of arguments)
        if (value.length <= 0) {
            if (brch[ibrch].type == "select") {
                if (iarg == 0) {
                    alert(name+" should not be blank (type)");
                    return;
                }
            } else if (brch[ibrch].type == "udparg" || brch[ibrch].type == "udprim") {
                if (iarg == 0) {
                    alert(name+" should not be blank (primtype)");
                    return;
                } else if (iarg == 2 || iarg == 4 || iarg == 6 || iarg == 8) {
                    if (prev_value.length > 0) {
                        alert(name+" should not be blank (follows non-blank)");
                        return;
                    }
                }
            } else {
                alert(name+" should not be blank");
                return;
            }
        }

        // prepend a dollar sign if a string argument (ie, the name starts with a dollar sign)
        if (name.charAt(0) != "$" || value.length <= 0) {
            output =       value;
        } else {
            output = "$" + value;
        }

        // check for bars
        if (value.indexOf("|") >= 0) {
            alert(name+" should not contain a bar");
            return;
        }

        // check if there were any changes
        if (output != brch[ibrch].args[iarg]) {
            if (wv.curMode == 3) {
                postMessage("Changing "+name+" from '"+brch[ibrch].args[iarg]+
                            "' to '"+output+
                            "' ====> Re-build is needed <====");
            }

            brch[ibrch].args[iarg] = output;
            nchange++;
        }

        mesg = mesg + output + "|";

        // save previous value
        prev_value = value;
    }

    // put rest of bars at end of mesg
    for (iarg = wv.numArgs; iarg < 9; iarg++) {
        mesg = mesg + "|";
    }

    if (wv.curMode == 3) {
        mesg = mesg + "||";
    }

    // if adding a new Branch, send the "newBrch" message
    if (wv.curMode == 2) {
        postMessage("Branch (type="+brch[ibrch].type+") has been added"+
                "  ====> Re-build is needed <====");

        browserToServer(mesg);

    // if there are changes, sent the "setBrch" message
    } else if (nchange > 0) {
        browserToServer(mesg);
    }

    // update any necessary attributes
    if (wv.curMode != 2) {
        for (var iattr = 0; iattr < brch[ibrch].attrs.length; iattr++) {
            if        (iattr == 0) {
                var name  = document.getElementById("attrName1").firstChild["data"];
                var value = editBrchForm.            attrValu1.value.replace(/\s/g, "");
            } else if (iattr == 1) {
                var name  = document.getElementById("attrName2").firstChild["data"];
                var value = editBrchForm.            attrValu2.value.replace(/\s/g, "");
            } else if (iattr == 2) {
                var name  = document.getElementById("attrName3").firstChild["data"];
                var value = editBrchForm.            attrValu3.value.replace(/\s/g, "");
            } else if (iattr == 3) {
                var name  = document.getElementById("attrName4").firstChild["data"];
                var value = editBrchForm.            attrValu4.value.replace(/\s/g, "");
            } else if (iattr == 4) {
                var name  = document.getElementById("attrName5").firstChild["data"];
                var value = editBrchForm.            attrValu5.value.replace(/\s/g, "");
            } else if (iattr == 5) {
                var name  = document.getElementById("attrName6").firstChild["data"];
                var value = editBrchForm.            attrValu6.value.replace(/\s/g, "");
            } else if (iattr == 6) {
                var name  = document.getElementById("attrName7").firstChild["data"];
                var value = editBrchForm.            attrValu7.value.replace(/\s/g, "");
            } else if (iattr == 7) {
                var name  = document.getElementById("attrName8").firstChild["data"];
                var value = editBrchForm.            attrValu8.value.replace(/\s/g, "");
            } else if (iattr == 8) {
                var name  = document.getElementById("attrName9").firstChild["data"];
                var value = editBrchForm.            attrValu9.value.replace(/\s/g, "");
            }

            // check for bars
            if (value.indexOf("|") >= 0) {
                alert(name+" should not contain a bar");
                return;
            }

            // check if there were any changes
            if (value != brch[ibrch].attrs[iattr][2]) {
                if        (brch[ibrch].attrs[iattr][1] == "(attr)") {
                    if (value.length > 0) {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|1|"+value+"|");
                        postMessage("Changing attribute '"+name+"' to "+
                                    value+" ====> Re-build is needed <====");
                    } else {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|1|<DeLeTe>|");
                        postMessage("Attribute '"+name+
                                    "' being deleted. ====> Re-build is needed <====");
                    }
                } else if (brch[ibrch].attrs[iattr][1] == "(csys)") {
                    if (value.length > 0) {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|2|"+value+"|");
                        postMessage("Changing csystem '"+name+"' to "+
                                    value+" ====> Re-build is needed <====");
                    } else {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|2|<DeLeTe>|");
                        postMessage("Csystem '"+name+
                                    "' being deleted. ====> Re-build is needed <====");
                    }
                } else {
                    alert("ERROR:: type="+brch[ibrch].attrs[iattr][1]);
                }

                brch[ibrch].attrs[iattr][2] = value;
                nchange++;
            }
        }
    }

    if (nchange > 0 || wv.curMode == 2) {
        wv.nchanges += nchange;

        // get an updated Branch list (so that added Branch is listed)
        browserToServer("getBrchs|");

        // update the UI
        if (wv.curMode == 3) {
            var id     = wv.menuEvent["target"].id;
            var myElem = document.getElementById(id);
            myElem.className = "fakelinkoff";
        }

        activateBuildButton();
    } else {
        alert("no changes were made");
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if not adding skbeg/skend, return to the WebViewer
    if (wv.curMode != 2 || brch[ibrch].type != "skbeg") {
        changeMode(0);

    //otherwise enter the sketcher
    } else {
        sket.ibrch = ibrch + 1;

        // send the message to the server
        browserToServer("loadSketch|"+sket.ibrch+"|");

        // inactivate buttons until build is done
        changeMode(-1);
    }
}


//
// callback when "Cancel" is pressed in editBrchForm (called by ESP.html)
//
function editBrchCancel() {
    // alert("in editBrchCancel()");

    // if we are in process of adding a Branch, remove it now
    if (wv.curMode == 2) {
        brch.splice(brch.length-1, 1);
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
}


//
// callback when Body name is pressed in Tree
//
function showBodyAttrs(e) {
    // alert("in showBodyAttrs()");

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    var myElem = document.getElementById("node"+inode+"col2");
    var bodyName = myElem.firstChild.nodeValue;
    while (bodyName.charCodeAt(0) == 160) {
        bodyName = bodyName.slice(1);
    }

    var mesg  = bodyName+":";
    try {
        var attrs = wv.sgData[bodyName];
        for (var i = 0; i < attrs.length; i+=2) {
            mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
        }
    } catch (x) {
    }
    postMessage(mesg);
}


//
// callback when any mouse is pressed in canvas (when wv.curMode==0)
//
function getMouseDown0(e)
{
    if (!e) var e = event;

    wv.startX   =  e.clientX - wv.offLeft             - 1;
    wv.startY   = -e.clientY + wv.offTop  + wv.height + 1;

    wv.dragging = true;
    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when the mouse moves in canvas (when wv.curMode==0)
//
function getMouseMove0(e)
{
    if (!e) var e = event;

    wv.cursorX  =  e.clientX - wv.offLeft             - 1;
    wv.cursorY  = -e.clientY + wv.offTop  + wv.height + 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when the mouse is released in canvas (when wv.curMode==0)
//
function getMouseUp0(e)
{
    wv.dragging = false;
}


//
// callback when the mouse wheel is rolled in canvas (when wv.curMode==0)
//
function getMouseRoll0(e)
{
    if (e) {

        // zoom in
        if        (e.deltaY > 0) {
            wv.mvMatrix.scale(1.1, 1.1, 1.1);
            wv.scale *= 1.1;
            wv.sceneUpd = 1;

        // zoom out
        } else if (e.deltaY < 0) {
            wv.mvMatrix.scale(0.9, 0.9, 0.9);
            wv.scale *= 0.9;
            wv.sceneUpd = 1;
        }
    }
}


//
// callback when the mouse leaves the canvas (when wv.curMode==0)
//
function mouseLeftCanvas0(e)
{
    if (wv.dragging) {
        wv.dragging = false;
    }
}


//
// callback when any mouse is pressed in Sketcher (when wv.curMode==7)
//
function getMouseDown7(e)
{
    if (!e) var e = event;

    wv.startX   = e.clientX - wv.offLeft7 - 1;
    wv.startY   = e.clientY - wv.offTop7  - 1;

    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    // set default key presses that are bound to mouse clicks in Sketcher
    if        (sket.mode == 1) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 108;   // "l"

        getKeyPress(evnt);
    } else if (sket.mode == 2) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 108;   // "l"

        getKeyPress(evnt);
    } else if (sket.mode == 3) {
        wv.cursorX = wv.startX;
        wv.cursorY = wv.startY;

        sket.movingPoint = getClosestSketchPoint();
    } else if (sket.mode == 4) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 124;   // "w"

        getKeyPress(evnt);
    } else if (sket.mode == 5) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 100;   // "d"

        getKeyPress(evnt);
    }
}


//
// callback when the mouse moves in Sketcher (when wv.curMode==7)
//
function getMouseMove7(e)
{
    if (!e) var e = event;

    wv.cursorX  = e.clientX - wv.offLeft7 - 1;
    wv.cursorY  = e.clientY - wv.offTop7  - 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    // if sket.mode==1, then the "proposed" new Segment
    //    is taken care of in drawSketch
    if (sket.mode == 1) {
        drawSketch();

    // if sket.mode==3, then move the point at (wv.startX, wv.startY)
    //    to the current location
    } else if (sket.mode == 3) {
        if (sket.movingPoint >= 0) {
            sket.pnt.x[sket.movingPoint] = wv.cursorX;
            sket.pnt.y[sket.movingPoint] = wv.cursorY;

            for (var iseg = 0; iseg < sket.seg.type.length; iseg++) {
                updateSegData(iseg);
            }

            drawSketch();
        }

    // if sket.mode==4 or 5, then the "proposed" new W or D constraint
    //    is taken care of in drawSketch
    } else if (sket.mode == 4 || sket.mode == 5) {
        drawSketch();

    // if sket.mode==2, then we need to use the cursor location to update
    //    the dip for the last Segment
    } else if (sket.mode == 2) {
        var iseg = sket.seg.type.length - 1;

        var xa = sket.pnt.x[sket.seg.ibeg[iseg]];
        var ya = sket.pnt.y[sket.seg.ibeg[iseg]];
        var xb = sket.pnt.x[sket.seg.iend[iseg]];
        var yb = sket.pnt.y[sket.seg.iend[iseg]];

        // put the cursor on the arc
        var xe  = wv.cursorX;
        var ye  = wv.cursorY;

        var D = (ya - ye) * (xb - xe) - (yb - ye) * (xa - xe);
        if (Math.abs(D) > 1e-6) {
            var s = ((xb - xa) * (xb - xe) - (ya - yb) * (yb - ye)) / D;

            var xc = (xa + xe + s * (ya - ye)) / 2;
            var yc = (ya + ye + s * (xe - xa)) / 2;

            var R = Math.sqrt((xc-xa) * (xc-xa) + (yc-ya) * (yc-ya));
            var L = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));

            if ((xe-xa)*(yb-ya) > (ye-ya)*(xb-xa)) {
                if ((xb-xa)*(yc-ya) > (yb-ya)*(xc-xa)) {
                    sket.seg.dip[iseg] = -R + Math.sqrt(R*R - L*L/4);
                } else {
                    sket.seg.dip[iseg] = -R - Math.sqrt(R*R - L*L/4);
                }
            } else {
                if ((xb-xa)*(yc-ya) > (yb-ya)*(xc-xa)) {
                    sket.seg.dip[iseg] = +R + Math.sqrt(R*R - L*L/4);
                } else {
                    sket.seg.dip[iseg] = +R - Math.sqrt(R*R - L*L/4);
                }
            }

            updateSegData(iseg);
        } else {
            sket.seg.dip[ iseg] = 0;
            sket.seg.xm[  iseg] = (xa + xb) / 2;
            sket.seg.ym[  iseg] = (ya + yb) / 2;
            sket.seg.xc[  iseg] = 0;
            sket.seg.yc[  iseg] = 0;
            sket.seg.rad[ iseg] = 0;
            sket.seg.tbeg[iseg] = 0;
            sket.seg.tend[iseg] = 0;
        }

        drawSketch();
    }
}


//
// callback when the mouse is released in Sketcher (when wv.curMode==7)
//
function getMouseUp7(e)
{
    sket.movingPoint = -1;
}


//
// callback when a key is pressed
//
function getKeyPress(e)
{
//    if (!e) var e = event;

    // if <esc> was pressed, return to base mode (if not in Sketcher)
    if (e.charCode == 0 && e.keyCode == 27 && wv.curMode != 7) {
        activateBuildButton();
        changeMode(0);
    }

    // if in canvas or Sketcher, record info about keypress
    if (wv.curMode == 0 || wv.curMode == 7) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

                        wv.modifier  = 0;
        if (e.shiftKey) wv.modifier |= 1;
        if (e.altKey  ) wv.modifier |= 2;
        if (e.ctrlKey ) wv.modifier |= 4;

        // if in the Sketcher, process the key press now
        if (wv.curMode == 7) {
            sketchKeyPress();
        }

    // if addBrchForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 1) {
        if (e.keyCode == 13) {
            addBrchOk();
            return false;
        }

    // if editBrchForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 2 || wv.curMode == 3) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

        if (wv.keyCode == 13) {
            editBrchOk();
            return false;
        } else {
            return true;
        }

    // if editPmtrForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 4 || wv.curMode == 5) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

        if (wv.keyCode == 13) {
            editPmtrOk();
            return false;
        } else {
            return true;
        }
    }

    return true;
}


//
// callback when an arrow... or shift is pressed (needed for Chrome)
//
function getKeyDown(e) {
    if (e.charCode == 0) {
        if (e.keyCode == 33 ||          // PgUp
            e.keyCode == 34 ||          // PgDn
            e.keyCode == 36 ||          // Home
            e.keyCode == 37 ||          // Left
            e.keyCode == 38 ||          // Up
            e.keyCode == 39 ||          // Right
            e.keyCode == 40   ) {       // Down
            wv.keyCode  = e.keyCode;
            wv.keyPress = 0;
        } else if (e.keyCode == 16) {   // Shift
            wv.modifier = 1;
        }
    }
}


//
// callback when a shift is released (needed for Chrome)
//
function getKeyUp(e) {

    if (e.charCode == 0 && e.keyCode == 16) {
        wv.modifier = 0;
    }
}


//
// callback when the mouse is pressed in key window
//
function setKeyLimits(e)
{

    // get new limits
    var templo = prompt("Enter new lower limit", wv.loLimit);
    if (isNaN(templo)) {
        alert("Lower limit must be a number");
        return;
    }

    var tempup = prompt("Enter new upper limit", wv.upLimit);
    if (isNaN(tempup)) {
        alert("Upper limit must be a number");
        return;
    }

    if (Number(tempup) <= Number(templo)) {
        alert("Upper limit must be greater than lower limit");
        return;
    }

    if (templo != wv.loLimit || tempup != wv.uplimit) {
        wv.loLimit = templo;
        wv.upLimit = tempup;

        // send the limits back to the server
        browserToServer("setLims|"+wv.plotType+"|"+wv.loLimit+"|"+wv.upLimit+"|");
    }
}


//
// callback to toggle Viz property
//
function toggleViz(e) {
    // alert("in toggleViz(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Viz property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
            myTree.prop(inode, 1, "on");
        } else {
            myTree.prop(inode, 1, "off");
        }

    //  toggle the Viz property (on all Faces/Edges in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col3");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 1, "on");
        } else {
            myTree.prop(inode, 1, "off");
        }
    }

    myTree.update();
}


//
// callback to toggle Grd property
//
function toggleGrd(e) {
    // alert("in toggleGrd(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Grd property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
            myTree.prop(inode, 2, "on");
        } else {
            myTree.prop(inode, 2, "off");
        }

    // toggle the Grd property (on all Faces/Edges in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col4");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 2, "on");
        } else {
            myTree.prop(inode, 2, "off");
        }
    }

    myTree.update();
}


//
// callback to toggle Trn property
//
function toggleTrn(e) {
    // alert("in toggleTrn(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Trn property (on a Face)
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
            myTree.prop(inode, 3, "on");
        } else {
            myTree.prop(inode, 3, "off");
        }

    // toggle the Trn property (on all Faces in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col5");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 3, "on");
            myElem.setAttribute("class", "fakelinkon");
            myElem.title = "Toggle Trn off";
        } else {
            myTree.prop(inode, 3, "off");
            myElem.setAttribute("class", "fakelinkoff");
            myElem.title = "Toggle Trn on";
        }
    }

    myTree.update();
}


//
// callback to toggle Ori property
//
function toggleOri(e) {
    // alert("in toggleOri(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Ori property (on an Edge)
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
            myTree.prop(inode, 3, "on");
        } else {
            myTree.prop(inode, 3, "off");
        }

    // toggle the Ori property (on all Edges in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col5");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 3, "on");
            myElem.setAttribute("class", "fakelinkon");
            myElem.title = "Toggle Ori off";
        } else {
            myTree.prop(inode, 3, "off");
            myElem.setAttribute("class", "fakelinkoff");
            myElem.title = "Toggle Ori on";
        }
    }

    myTree.update();
}


//
// callback when "CancelStepThru" is pressed in Tree
//
function cancelStepThru() {
    // alert("in cancelStepThru()");

    wv.curStep = 0;
    var button = document.getElementById("stepThruBtn");
    button["innerHTML"] = "StepThru";

    browserToServer("nextStep|0|");
}


//
// callback when "DisplayType" is pressed
//
function modifyDisplayType(e) {
    // alert("in modifyDisplayType(e="+e+")");

    var ptype = prompt("Enter display type:\n"+
                      "   0  monochrome\n"+
                      "   1  normalized U parameter\n"+
                      "   2  normalized V parameter\n"+
                      "   3  minimum curvature\n"+
                      "   4  maximum curvature\n"+
                      "   5  Gaussian curvature", "0");

    if (ptype === null) {
        return;
    } else if (isNaN(ptype)) {
        alert("Illegal number format, so no change");
        return;
    } else {
        wv.plotType = Number(ptype);

        setKeyLimits(null);
    }
}


//
// callback when "DisplayFilter" is pressed
//
function modifyDisplayFilter(e) {
    // alert("in modifyDisplayFilter(e="+e+")")

    var attrName = prompt("Enter Attribute name (or * for all)");
    if (attrName === null) {
        displayFilterOff();
    } else if (attrName.length <= 0) {
        displayFilterOff();
    } else {
        var attrValue = prompt("Enter Attribute value (or * for all)");
        if (attrValue === null) {
            displayFilterOff();
        } else if (attrValue.length <= 0) {
            displayFilterOff();
        } else {
            for (var gprim in wv.sceneGraph) {
                // make all gprims transparent and lWidth=2 for Edges
                wv.sceneGraph[gprim].attrs |= wv.plotAttrs.TRANSPARENT;
                if (wv.sceneGraph[gprim].lWidth == 5) {
                    wv.sceneGraph[gprim].lWidth = 2;
                }
                try {
                    // now make those that match the filter non-transparent
                    var attrs = wv.sgData[gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        if        (attrs[i].trim() == attrName && attrs[i+1].trim() == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            if (wv.sceneGraph[gprim].lWidth == 2) {
                                wv.sceneGraph[gprim].lWidth = 5;
                            }
                        } else if ("*"             == attrName && attrs[i+1].trim() == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            if (wv.sceneGraph[gprim].lWidth == 2) {
                                wv.sceneGraph[gprim].lWidth = 5;
                            }
                        } else if (attrs[i].trim() == attrName && "*"               == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            if (wv.sceneGraph[gprim].lWidth == 2) {
                                wv.sceneGraph[gprim].lWidth = 5;
                            }
                        } else if ("*"             == attrName && "*"               == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            if (wv.sceneGraph[gprim].lWidth == 2) {
                                wv.sceneGraph[gprim].lWidth = 5;
                            }
                        }
                    }
                } catch (x) {
                }
            }
            postMessage("Display filtered to \""+attrName+"\" \""+attrValue+"\"");
        }
    }

    myTree.update();
    wv.sceneUpd = 1;
}



//
// turn dsiplay filtering off
//
function displayFilterOff() {
    // alert("in displayFilterOff()");

    for (var gprim in wv.sceneGraph) {
        wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
        if (wv.sceneGraph[gprim].lWidth == 5) {
            wv.sceneGraph[gprim].lWidth = 2;
        }
    }

    postMessage("Display filtering off");
}


//
// Initialize a Sketch
//
function initializeSketch()
{
    // alert("initializeSketch()");

    if (sket) {
        sket.basePoint = undefined;
        sket.skbegX    = undefined;
        sket.skbegY    = undefined;
        sket.skbegZ    = undefined;
        sket.relative  =  0;
        sket.suggest   = undefined;
        sket.halo      =  5;      // pixels for determining "closeness" while drawing
        wv.keyPress    = -1;      // last key pressed
        wv.keyCode     = -1;
        wv.button      = -1;      // button pressed
        wv.modifier    =  0;      // modifier (shift,alt,ctrl) bitflag

        // set up Sketcher offsets
        var sketcher = document.getElementById("sketcher");
        wv.offTop7   = sketcher.offsetTop;
        wv.offLeft7  = sketcher.offsetLeft;

        // initialize the arrays
        sket.pnt        = {};
        sket.pnt.x      = [];
        sket.pnt.y      = [];
        sket.pnt.lbl    = [];

        sket.seg        = {};
        sket.seg.type   = [];
        sket.seg.ibeg   = [];
        sket.seg.iend   = [];
        sket.seg.xm     = [];
        sket.seg.ym     = [];
        sket.seg.lbl    = [];

        sket.seg.dip    = [];
        sket.seg.xc     = [];
        sket.seg.yc     = [];
        sket.seg.rad    = [];
        sket.seg.tbeg   = [];
        sket.seg.tend   = [];

        sket.var        = {};
        sket.var.name   = [];
        sket.var.value  = [];

        sket.con        = {};
        sket.con.type   = [];
        sket.con.index1 = [];
        sket.con.index2 = [];
        sket.con.value  = [];

        // mapping between "physical" and "screen" coordinates
        //   x_physical = sket.scale * (x_screen - sket.xorig)
        //   y_physical = sket.scale * (sket.yorig - y_screen)

        // commands that set sket.xorig: "x"
        // commands that set sket.yorig: "y"
        // commands that set sket.scale: "l", "r", "w", or "d"
        //                               "x" if another "x" exists
        //                               "y" if another "y" exists

        // reset scaling info
        sket.scale = undefined;
        sket.xorig = undefined;
        sket.yorig = undefined;

        // initialize undo
        sket.nundo   = 0;
        sket.undo    = new Array;
        for (var iundo = 0; iundo < sket.maxundo; iundo++) {
            sket.undo[iundo] = undefined;
        }

        // set mode to initializing
        sket.mode = 0;
    } else {
        alert("ERROR:: sket does not exist");
    }

    // draw the current Sketch
    drawSketch();
}


//
// load a Sketch from a JSON string sent from serveCSM
//
function loadSketch(begs, vars, cons, segs)
{
    // alert("loadSketch()");
    // alert("   begs="+begs);
    // alert("   vars="+vars);
    // alert("   cons="+cons);
    // alert("   segs="+segs);

    // break up the inputs
    var begsList = begs.split(";");
    var varsList = vars.split(";");
    var consList = cons.split(";");
    var segsList = segs.split(";");

    // get the number of points, segments, variables, and constraints
    var npnt = (varsList.length - 1) / 3;
    var ncon = (consList.length - 1) / 4;
    var nseg = (segsList.length - 1) / 3;

    // store the base expressions from the SKBEG statement
    sket.skbegX   = begsList[0];
    sket.skbegY   = begsList[2];
    sket.skbegZ   = begsList[4];
    sket.relative = begsList[6];
    sket.suggest  = undefined;

    // find the extrema of the coordinates
    var xmin = Number(begsList[1]);
    var xmax = Number(begsList[1]);
    var ymin = Number(begsList[3]);
    var ymax = Number(begsList[3]);

    for (var ipnt = 1; ipnt < npnt; ipnt++) {
        xmin = Math.min(xmin, Number(varsList[3*ipnt  ]));
        xmax = Math.max(xmax, Number(varsList[3*ipnt  ]));
        ymin = Math.min(ymin, Number(varsList[3*ipnt+1]));
        ymax = Math.max(ymax, Number(varsList[3*ipnt+1]));
    }

    // get the size of the canvas
    var canvas = document.getElementById("sketcher");
    var width  = canvas.clientWidth;
    var height = canvas.clientHeight;

    if (npnt > 1) {
        var test1 = (xmax - xmin) / width;
        var test2 = (ymax - ymin) / height;

        // set up sizing so that Sketch fills middle 50 percent of the window
        sket.scale = 2 * Math.max(test1, test2);
        sket.xorig = (width  - (xmin + xmax) / sket.scale) / 2;
        sket.yorig = (height + (ymin + ymax) / sket.scale) / 2;

    // not enough data to set scales, so guess something base upon the SKBEG statement
    } else {
        sket.scale = 0.01;
        sket.xorig = width  / 2 - begsList[1] / sket.scale;
        sket.yorig = height / 2 + begsList[3] / sket.scale;
    }

    // create the points
    sket.pnt.x   = [];
    sket.pnt.y   = [];
    sket.pnt.lbl = [];

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        sket.pnt.x[  ipnt] = Math.floor(sket.xorig + Number(varsList[3*ipnt  ]) / sket.scale);
        sket.pnt.y[  ipnt] = Math.floor(sket.yorig - Number(varsList[3*ipnt+1]) / sket.scale);
        sket.pnt.lbl[ipnt] = "";
    }

    // overwrite the first point by the begsList
    if (npnt > 0) {
        sket.pnt.x[0] = Math.floor(sket.xorig + Number(begsList[1]) / sket.scale);
        sket.pnt.y[0] = Math.floor(sket.yorig - Number(begsList[3]) / sket.scale);
    }

    // create the segments
    sket.seg.type = [];
    sket.seg.ibeg = [];
    sket.seg.iend = [];
    sket.seg.xm   = [];
    sket.seg.ym   = [];
    sket.seg.lbl  = [];

    sket.seg.dip  = [];
    sket.seg.xc   = [];
    sket.seg.yc   = [];
    sket.seg.rad  = [];
    sket.seg.tbeg = [];
    sket.seg.tend = [];

    for (var iseg = 0; iseg < nseg; iseg++) {
        sket.seg.type[iseg] =        segsList[3*iseg  ];
        sket.seg.ibeg[iseg] = Number(segsList[3*iseg+1]) - 1;
        sket.seg.iend[iseg] = Number(segsList[3*iseg+2]) - 1;
        sket.seg.lbl[ iseg] = "";

        sket.seg.dip[ iseg] = 0;
        sket.seg.xc[  iseg] = 0;
        sket.seg.yc[  iseg] = 0;
        sket.seg.rad[ iseg] = 0;
        sket.seg.tbeg[iseg] = 0;
        sket.seg.tend[iseg] = 0;

        // update the segment if an arc
        if (sket.seg.type[iseg] == "C") {
            if (iseg < npnt-1) {
                sket.seg.dip[iseg] = Number(varsList[3*iseg+5]) / sket.scale;
            } else {
                sket.seg.dip[iseg] = Number(varsList[       2]) / sket.scale;
            }
        }

        updateSegData(iseg);
    }

    // create the variables
    sket.var.name  = [];
    sket.var.value = [];

    var nvar = 0;
    for (var ivar = 0; ivar < npnt; ivar++) {
        sket.var.name[ nvar] = sprintf("::x[%]", ivar+1);
        sket.var.value[nvar] = varsList[3*ivar  ];
        nvar++;

        sket.var.name[ nvar] = sprintf("::y[%]", ivar+1);
        sket.var.value[nvar] = varsList[3*ivar+1];
        nvar++;

        if (varsList.length > 0 && Number(varsList[3*ivar+2]) != 0) {
            sket.var.name[ nvar] = sprintf("::d[%]", ivar+1);
            sket.var.value[nvar] = varsList[3*ivar+2];
            nvar++;
        }
    }

    // create the constraints
    sket.con.type   = [];
    sket.con.index1 = [];
    sket.con.index2 = [];
    sket.con.value  = [];

    for (var icon = 0; icon < ncon; icon++) {
        sket.con.type[  icon] =        consList[4*icon  ];
        sket.con.index1[icon] = Number(consList[4*icon+1]) - 1;
        sket.con.index2[icon] = Number(consList[4*icon+2]) - 1;
        sket.con.value[ icon] =        consList[4*icon+3];

        if (sket.con.type[icon] == "W" || sket.con.type[icon] == "D") {
            // not associated with Point or Segment
        } else if (sket.con.index2[icon] < 0) {
            sket.con.index2[icon] += 1;
            sket.pnt.lbl[sket.con.index1[icon]] += sket.con.type[icon];
        } else {
            sket.seg.lbl[sket.con.index1[icon]] += sket.con.type[icon];
        }
    }

    // if there were no constraints, create X and Y constraints at first Point
    if (ncon == 0) {
        sket.pnt.x[  0] = Math.floor(sket.xorig + Number(begsList[1]) / sket.scale);
        sket.pnt.y[  0] = Math.floor(sket.yorig - Number(begsList[3]) / sket.scale);
        sket.pnt.lbl[0] = "XY";

        sket.var.name[ 0] = "::x[1]";
        sket.var.value[0] = sket.pnt.x[0];

        sket.var.name[ 1] = "::y[1]";
        sket.var.value[1] = sket.pnt.y[0];

        if (sket.relative == 0) {
            sket.con.type[  0] = "X";
            sket.con.index1[0] = 0;
            sket.con.index2[0] = -1;
            sket.con.value[ 0] = sket.skbegX;

            sket.con.type[  1] = "Y";
            sket.con.index1[1] = 0;
            sket.con.index2[1] = -1;
            sket.con.value[ 1] = sket.skbegY;
        } else {
            sket.con.type[  0] = "X";
            sket.con.index1[0] = 0;
            sket.con.index2[0] = -1;
            sket.con.value[ 0] = 0;

            sket.con.type[  1] = "Y";
            sket.con.index1[1] = 0;
            sket.con.index2[1] = -1;
            sket.con.value[ 1] = 0;
        }
    }

    // set the mode depending on status of Sketch
    if (sket.seg.type.length == 0) {
        sket.mode = 1;
    } else if (sket.seg.iend[sket.seg.iend.length-1] != 0) {
        sket.mode = 1;
    } else {
        sket.mode = 3;
    }

    // draw the current Sketch
    drawSketch();
}


//
// process a key press in the Sketcher
//
function sketchKeyPress()
{
    // alert("sketchKeyPress");

    var canvas  = document.getElementById("sketcher");
    var context = canvas.getContext("2d");

    var npnt = sket.pnt.x.length;
    var nseg = sket.seg.type.length;
    var nvar = sket.var.name.length;
    var ncon = sket.con.type.length;

    var myKeyPress = String.fromCharCode(wv.keyPress);

    // 'ctrl-h' --- home
    if ((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode == 0) ||
        (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode == 8)   ) {
        cmdHome();
        return;

    // '<Home>' -- home
    } else if (wv.keyPress == 0 && wv.keyCode == 36) {
        cmdHome();
        return;

    // 'ctrl-i' - zoom in (same as <PgUp> without shift)
    } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
        cmdIn();
        return;

    // '<PgDn>' -- zoom in
    } else if (wv.keyPress == 0 && wv.keyCode == 33) {
        cmdIn();
        return;

    // '+' -- zoom in
    } else if (wv.keyPress == 43 && wv.modifier == 1) {
        cmdIn();
        return;

    // 'ctrl-o' - zoom out (same as <PgDn> without shift)
    } else if ((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) {
        cmdOut();
        return;

    // '<PgUp>' -- zoom out
    } else if (wv.keyPress == 0 && wv.keyCode == 34) {
        cmdOut();
        return;

    // '-' -- zoomout
    } else if (wv.keyPress == 45 && wv.modifier == 0) {
        cmdOut();
        return;

    // 'ctrl-r' - translate right
    } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
        cmdRite();
        return;

    // '<Right>' -- translate right
    } else if (wv.keyPress == 0 && wv.keyCode == 39) {
        cmdRite();
        return;

    // 'ctrl-l' - translate left
    } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
        cmdLeft();
        return;

    // '<Left>' -- translate left
    } else if (wv.keyPress == 0 && wv.keyCode == 37) {
        cmdLeft();
        return;

    // 'ctrl-t' - translate up
    } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
        cmdTop();
        return;

    // '<Up>' -- translate up
    } else if (wv.keyPress == 0 && wv.keyCode == 38) {
        cmdTop();
        return;

    // 'ctrl-b' - translate down
    } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
        cmdBotm();
        return;

    // '<Down>' -- translate down
    } else if (wv.keyPress == 0 && wv.keyCode == 40) {
        cmdBotm();
        return;

    // '@' - coordinates of cursor
    } else if (myKeyPress == "@") {
        var xx = sket.scale * (wv.cursorX - sket.xorig);
        var yy = sket.scale * (sket.yorig - wv.cursorY);
        postMessage("Cursor is at ("+xx+","+yy+")");

    // '*' - rescale based upon best info
    } else if (myKeyPress == '*') {
        var xold = sket.scale * (wv.cursorX - sket.xorig);
        var yold = sket.scale * (sket.yorig - wv.cursorY);

        var xnew = prompt("Enter approximate X location at cursor: ", xold);
        if (xnew !== null) {
            var ynew = prompt("Enter approximate Y location at cursor: ", yold);
            if (ynew !== null) {
                postMessage("resetting scale based upon cursor location");

                var Lphys = Math.sqrt(xnew*xnew + ynew*ynew);
                var Lscr  = Math.sqrt(Math.pow(wv.cursorX-sket.xorig,2)
                                     +Math.pow(wv.cursorY-sket.yorig,2));
                sket.scale = Lphys / Lscr;
            }
        }

    // sket.mode 1 options ('l', 'c', 's', 'b', 'z', or 'o' to leave sketch open)
    } else if (sket.mode == 1) {
        saveSketchUndo();

        // determine the coordinates (and constraints) for this placement
        var xnew = wv.cursorX;
        var ynew = wv.cursorY;
        var cnew = "";

        for (var ipnt = npnt-1; ipnt >= 0; ipnt--) {
            if (Math.abs(wv.cursorX-sket.pnt.x[ipnt]) < sket.halo) {
                xnew = sket.pnt.x[ipnt];
                if (ipnt == npnt-1) {
                    cnew = "V";
                }
                break;
            }
        }
        for (ipnt = npnt-1; ipnt >= 0; ipnt--) {
            if (Math.abs(wv.cursorY-sket.pnt.y[ipnt]) < sket.halo) {
                ynew = sket.pnt.y[ipnt];
                if (ipnt == npnt-1) {
                    cnew = "H";
                }
                break;
            }
        }

        // "l" --- add a line Segment
        if (myKeyPress == "l" || myKeyPress == "L") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of last point, reject it
            if (Math.abs(sket.pnt.x[npnt-1]-xnew) < sket.halo &&
                Math.abs(sket.pnt.y[npnt-1]-ynew) < sket.halo   ) {
                postMessage("cannot create zero-length line segment");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sket.pnt.x[0]-xnew) > sket.halo ||
                Math.abs(sket.pnt.y[0]-ynew) > sket.halo   ) {
                iend = npnt;
                sket.pnt.x[  iend] = xnew;
                sket.pnt.y[  iend] = ynew;
                sket.pnt.lbl[iend] = "";

                sket.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sket.var.value[nvar  ] = sprintf("%",    xnew  );

                sket.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sket.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sket.seg.type[nseg] = "L";
            sket.seg.ibeg[nseg] = ibeg;
            sket.seg.iend[nseg] = iend;
            sket.seg.xm[  nseg] = (sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 2;
            sket.seg.ym[  nseg] = (sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 2;
            sket.seg.lbl[ nseg] = "";

            sket.seg.dip[ nseg] = 0;
            sket.seg.xc[  nseg] = 0;
            sket.seg.yc[  nseg] = 0;
            sket.seg.rad[ nseg] = 0;
            sket.seg.tbeg[nseg] = 0;
            sket.seg.tend[nseg] = 0;

            // add a new "horizontal" or "vertical" constraint
            if        (cnew == "H") {
                sket.con.type[  ncon] = "H";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[nseg] += "H";
            } else if (cnew == "V") {
                sket.con.type[  ncon] = "V";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[nseg] += "V";
            }

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sket.mode = 3;
            }

        // "c" --- add a circular arc
        } else if (myKeyPress == "c" || myKeyPress == "C") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of last point, reject it
            if (Math.abs(sket.pnt.x[npnt-1]-xnew) < sket.halo &&
                Math.abs(sket.pnt.y[npnt-1]-ynew) < sket.halo   ) {
                postMessage("cannot create zero-length circular arc");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sket.pnt.x[0]-xnew) > sket.halo ||
                Math.abs(sket.pnt.y[0]-ynew) > sket.halo) {
                iend          = npnt;
                sket.pnt.x[  iend] = xnew;
                sket.pnt.y[  iend] = ynew;
                sket.pnt.lbl[iend] = "";

                sket.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sket.var.value[nvar  ] = sprintf("%",    xnew  );

                sket.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sket.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sket.seg.type[nseg] = "C";
            sket.seg.ibeg[nseg] = ibeg;
            sket.seg.iend[nseg] = iend;
            sket.seg.xm[  nseg] = (sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 2;
            sket.seg.ym[  nseg] = (sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 2;
            sket.seg.lbl[ nseg] = "";

            sket.seg.dip[ nseg] = 0;
            sket.seg.xc[  nseg] = 0;
            sket.seg.yc[  nseg] = 0;
            sket.seg.rad[ nseg] = 0;
            sket.seg.tbeg[nseg] = 0;
            sket.seg.tend[nseg] = 0;

            nvar = sket.var.name.length;

            sket.var.name[ nvar] = sprintf("::d[%]", sket.seg.type.length);
            sket.var.value[nvar] = "0";

            sket.mode = 2;

        // "s" --- add spline Point
        } else if (myKeyPress == "s" || myKeyPress == "S") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of last point, reject it
            if (Math.abs(sket.pnt.x[npnt-1]-xnew) < sket.halo &&
                Math.abs(sket.pnt.y[npnt-1]-ynew) < sket.halo   ) {
                postMessage("cannot create zero-length spline segment");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sket.pnt.x[0]-xnew) > sket.halo ||
                Math.abs(sket.pnt.y[0]-ynew) > sket.halo   ) {
                iend = npnt;
                sket.pnt.x[  iend] = xnew;
                sket.pnt.y[  iend] = ynew;
                sket.pnt.lbl[iend] = "";

                sket.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sket.var.value[nvar  ] = sprintf("%",    xnew  );

                sket.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sket.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sket.seg.type[nseg] = "S";
            sket.seg.ibeg[nseg] = ibeg;
            sket.seg.iend[nseg] = iend;
            sket.seg.xm[  nseg] = (sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 2;
            sket.seg.ym[  nseg] = (sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 2;
            sket.seg.lbl[ nseg] = "";

            sket.seg.dip[ nseg] = 0;
            sket.seg.xc[  nseg] = 0;
            sket.seg.yc[  nseg] = 0;
            sket.seg.rad[ nseg] = 0;
            sket.seg.tbeg[nseg] = 0;
            sket.seg.tend[nseg] = 0;

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sket.mode = 3;
            }

        // "b" --- add Bezier Point
        } else if (myKeyPress == "b" || myKeyPress == "B") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of last point, reject it
            if (Math.abs(sket.pnt.x[npnt-1]-xnew) < sket.halo &&
                Math.abs(sket.pnt.y[npnt-1]-ynew) < sket.halo   ) {
                postMessage("cannot create zero-length bezier segment");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sket.pnt.x[0]-xnew) > sket.halo ||
                Math.abs(sket.pnt.y[0]-ynew) > sket.halo   ) {
                iend = npnt;
                sket.pnt.x[  iend] = xnew;
                sket.pnt.y[  iend] = ynew;
                sket.pnt.lbl[iend] = "";

                sket.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sket.var.value[nvar  ] = sprintf("%",    xnew  );

                sket.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sket.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sket.seg.type[nseg] = "B";
            sket.seg.ibeg[nseg] = ibeg;
            sket.seg.iend[nseg] = iend;
            sket.seg.xm[  nseg] = (sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 2;
            sket.seg.ym[  nseg] = (sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 2;
            sket.seg.lbl[ nseg] = "";

            sket.seg.dip[ nseg] = 0;
            sket.seg.xc[  nseg] = 0;
            sket.seg.yc[  nseg] = 0;
            sket.seg.rad[ nseg] = 0;
            sket.seg.tbeg[nseg] = 0;
            sket.seg.tend[nseg] = 0;

            // add a new "horizontal" or "vertical" constraint
            if        (cnew == "H") {
                sket.con.type[  ncon] = "H";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[nseg] += "H";
            } else if (cnew == "V") {
                sket.con.type[  ncon] = "V";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[nseg] += "V";
            }

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sket.mode = 3;
            }

        // "z" --- add a zero-length Segment (to break a spline or Bezier)
        } else if (myKeyPress == 'z') {
            var ibeg = npnt - 1;
            var iend = 0;

            // make sure this follows a Bezier or spline
            if (sket.seg.type[nseg-1] != "B" && sket.seg.type[nseg-1] != "S") {
                postMessage("zero-length segment can only come after spline or bezier");
                return;
            }

            // add a new Point and its variables
            sket.pnt.lbl[ibeg] = "Z";

            iend = npnt;
            sket.pnt.x[  iend] = sket.pnt.x[ibeg];
            sket.pnt.y[  iend] = sket.pnt.y[ibeg];
            sket.pnt.lbl[iend] = "Z";

            sket.var.name[ nvar  ] = sprintf("::x[%]", npnt);
            sket.var.value[nvar  ] = sprintf("%",    xnew  );

            sket.var.name[ nvar+1] = sprintf("::y[%]", npnt);
            sket.var.value[nvar+1] = sprintf("%",    ynew  );

            // add a new (linear) Segment
            sket.seg.type[nseg] = "L";
            sket.seg.ibeg[nseg] = ibeg;
            sket.seg.iend[nseg] = iend;
            sket.seg.xm[  nseg] = sket.pnt.x[ibeg];
            sket.seg.ym[  nseg] = sket.pnt.y[ibeg];
            sket.seg.lbl[ nseg] = "Z";

            sket.seg.dip[ nseg] = 0;
            sket.seg.xc[  nseg] = 0;
            sket.seg.yc[  nseg] = 0;
            sket.seg.rad[ nseg] = 0;
            sket.seg.tbeg[nseg] = 0;
            sket.seg.tend[nseg] = 0;

            // add the constraints
            sket.con.type[  ncon] = "Z";
            sket.con.index1[ncon] = ibeg;
            sket.con.index2[ncon] = -2;
            sket.con.value[ ncon] = 0;

            sket.con.type[  ncon+1] = "Z";
            sket.con.index1[ncon+1] = iend;
            sket.con.index2[ncon+1] = -3;
            sket.con.value[ ncon+1] = 0;

        // "o" --- end (open) Sketch
        } else if (myKeyPress == "o") {
            sket.mode = 3;

        } else {
            postMessage("Valid options are 'l', 'c', 's', 'b', 'z', and 'o'");
        }

    // sket.mode 2 option (any keyPress or mouseDown)
    } else if (sket.mode == 2) {

        // <any> --- set curvature
        if (sket.seg.iend[sket.seg.iend.length-1] == 0) {
            sket.mode = 3;
        } else {
            sket.mode = 1;
        }

    // sket.mode 3 options
    // if at a Point:          'x', 'y', 'p', 't', 'a', 'w', 'd' or '<' to delete
    // if at a cirarc Segment: 'r', 's', 'x', 'y',               or '<' to delete
    // if at          Segment: 'h', 'v', 'i', 'l',               or '<' to delete
    } else if (sket.mode == 3 || sket.mode == 6) {
        saveSketchUndo();

        // "x" --- set X coordinate of nearest Point or arc Segment
        if (myKeyPress == "x" || myKeyPress == "X") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            var jbest = getClosestSketchSegment();
            if (jbest >= 0 && sket.seg.type[jbest] != "C") {
                jbest = -1;
            }

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sket.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sket.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sket.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sket.seg.ym[jbest]);
            }

            var xvalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (dibest < djbest) {
                if (sket.pnt.lbl[ibest].indexOf("X") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sket.con.type[icon] == "X" && sket.con.index1[icon] == ibest) {
                            if (ibest == 0 && sket.relative == 1) {
                                alert("X cannot be changed at first point if in relative mode");
                            } else {
                                xvalue = prompt("Enter x value: ", sket.con.value[icon]);
                                if (xvalue !== null) {
                                    sket.con.value[icon] = xvalue;
                                }
                            }
                            break;
                        }
                    }
                } else {
                    if (sket.scale === undefined && sket.xorig === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sket.con.type[icon] == "X" && sket.con.index1[icon] != ibest) {
                                    sket.scale = (xvalue - sket.con.value[icon]) / (sket.pnt.x[ibest] - sket.pnt.x[sket.con.index1[icon]]);
                                    sket.xorig = sket.pnt.x[sket.con.index1[icon]] - sket.con.value[icon] / sket.scale;
                                    break;
                                }
                            }
                        }
                    } else if (sket.xorig === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            sket.xorig = sket.pnt.x[ibest] + xvalue / sket.scale;
                        }
                    } else if (sket.scale === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sket.con.type[icon] == "X" && sket.con.index1[icon] != ibest) {
                                    sket.scale = (xvalue - sket.con.value[icon]) / (sket.pnt.x[ibest] - sket.xorig);
                                    break;
                                }
                            }
                        }
                    } else {
                        xvalue = prompt("Enter x value: ", sket.scale*(sket.pnt.x[ibest]-sket.xorig));
                    }

                    if (xvalue !== null) {
                        sket.con.type[  ncon] = "X";
                        sket.con.index1[ncon] = ibest;
                        sket.con.index2[ncon] = -1;
                        sket.con.value[ ncon] = xvalue;

                        sket.pnt.lbl[ibest] += "X";
                    }
                }
            } else {
                if (sket.seg.lbl[jbest].indexOf("X") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sket.con.type[  icon] == "X" &&
                            sket.con.index1[icon] == sket.seg.ibeg[jbest] &&
                            sket.con.index2[icon] == sket.seg.iend[jbest]   ) {
                            xvalue = prompt("Enter x value: ", sket.con.value[icon]);
                            if (xvalue !== null) {
                                sket.con.value[icon] = xvalue;
                            }
                            break;
                        }
                    }
                } else {
                    xvalue = prompt("Enter x value: ");
                    if (xvalue !== null) {
                        sket.con.type[  ncon] = "X";
                        sket.con.index1[ncon] = sket.seg.ibeg[jbest];
                        sket.con.index2[ncon] = sket.seg.iend[jbest];
                        sket.con.value[ ncon] = xvalue;

                        sket.seg.lbl[jbest] += "X";
                    }
                }
            }

        // "y" --- set Y coordinate of nearest Point or arc Segment
        } else if (myKeyPress == "y" || myKeyPress == "Y") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            var jbest = getClosestSketchSegment();
            if (jbest >= 0 && sket.seg.type[jbest] != "C") {
                jbest = -1;
            }

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sket.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sket.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sket.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sket.seg.ym[jbest]);
            }

            var yvalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (dibest < djbest) {
                if (sket.pnt.lbl[ibest].indexOf("Y") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sket.con.type[icon] == "Y" && sket.con.index1[icon] == ibest) {
                            if (ibest == 0 && sket.relative == 1) {
                                alert("Y cannot be changed at first point if in relative mode");
                            } else {
                                yvalue = prompt("Enter y value: ", sket.con.value[icon]);
                                if (yvalue !== null) {
                                    sket.con.value[icon] = yvalue;
                                }
                            }
                            break;
                        }
                    }
                } else {
                    if (sket.scale === undefined && sket.yorig === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sket.con.type[icon] == "Y" && sket.con.index1[icon] != ibest) {
                                    sket.scale = (sket.con.value[icon] - yvalue) / (sket.pnt.y[ibest] - sket.pnt.y[sket.con.index1[icon]]);
                                    sket.yorig = sket.pnt.y[sket.con.index1[icon]] + sket.con.value[icon] / sket.scale;
                                    break;
                                }
                            }
                        }
                    } else if (sket.yorig === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            sket.yorig = sket.pnt.y[ibest] - yvalue / sket.scale;
                        }
                    } else if (sket.scale === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sket.con.type[icon] == "Y" && sket.con.index1[icon] != ibest) {
                                    sket.scale = (sket.con.value[icon] - yvalue) / (sket.pnt.y[ibest] - sket.yorig);
                                    break;
                                }
                            }
                        }
                    } else {
                        yvalue = prompt("Enter y value: ", sket.scale*(sket.yorig-sket.pnt.y[ibest]));
                    }

                    if (yvalue !== null) {
                        sket.con.type[  ncon] = "Y";
                        sket.con.index1[ncon] = ibest;
                        sket.con.index2[ncon] = -1;
                        sket.con.value[ ncon] = yvalue;

                        sket.pnt.lbl[ibest] += "Y";
                    }
                }
            } else {
                if (sket.seg.lbl[jbest].indexOf("Y") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sket.con.type[  icon] == "Y" &&
                            sket.con.index1[icon] == sket.seg.ibeg[jbest] &&
                            sket.con.index2[icon] == sket.seg.iend[jbest]   ) {
                            yvalue = prompt("Enter y value: ", sket.con.value[icon]);
                            if (xvalue !== null) {
                                sket.con.value[icon] = yvalue;
                            }
                            break;
                        }
                    }
                } else {
                    yvalue = prompt("Enter y value: ");
                    if (yvalue !== null) {
                        sket.con.type[  ncon] = "Y";
                        sket.con.index1[ncon] = sket.seg.ibeg[jbest];
                        sket.con.index2[ncon] = sket.seg.iend[jbest];
                        sket.con.value[ ncon] = yvalue;

                        sket.seg.lbl[jbest] += "Y";
                    }
                }
            }

        // "p" --- perpendicularity
        } else if (myKeyPress == "p" || myKeyPress == "P") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            if (ibest < 0) {
                alert("No Point found");
            } else if (sket.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("A") >= 0) {
                alert("'A' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else {
                sket.con.type[  ncon] = "P";
                sket.con.index1[ncon] = ibest;
                sket.con.index2[ncon] = -1;
                sket.con.value[ ncon] = 0;

                sket.pnt.lbl[ibest] += "P";
            }

        // "t" --- tangency
        } else if (myKeyPress == "t" || myKeyPress == "T") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            if (ibest < 0) {
                alert("No Point found");
            } else if (sket.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("A") >= 0) {
                alert("'A' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else {
                sket.con.type[  ncon] = "T";
                sket.con.index1[ncon] = ibest;
                sket.con.index2[ncon] = -1;
                sket.con.value[ ncon] = 0;

                sket.pnt.lbl[ibest] += "T";
            }

        // "a" --- angle (positive to left)
        } else if (myKeyPress == 'a' || myKeyPress == "A") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            var avalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (sket.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else if (sket.pnt.lbl[ibest].indexOf("A") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[icon] == "A" && sket.con.index1[icon] == ibest) {
                        avalue = prompt("Enter angle (deg): ", sket.con.value[icon]);
                        if (avalue !== null) {
                            sket.con.value[icon] = avalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = -1;
                var iend = -1;
                for (var iseg = 0; iseg < sket.seg.type.length; iseg++) {
                    if (sket.seg.iend[iseg] == ibest) {
                        ibeg = sket.seg.ibeg[iseg];
                    }
                    if (sket.seg.ibeg[iseg] == ibest) {
                        iend = sket.seg.iend[iseg];
                    }
                }

                if (ibeg < 0 || iend < 0) {
                    alert("Angle cannot be specified at endpoints");
                    return;
                }

                var ang1 = Math.atan2(sket.pnt.y[ibest]-sket.pnt.y[ibeg ],
                                      sket.pnt.x[ibest]-sket.pnt.x[ibeg ]) * 180 / Math.PI;
                var ang2 = Math.atan2(sket.pnt.y[iend ]-sket.pnt.y[ibest],
                                      sket.pnt.x[iend ]-sket.pnt.x[ibest]) * 180 / Math.PI;
                avalue = ang1 - ang2;
                while (avalue < -180) {
                    avalue += 360;
                }
                while (avalue > 180) {
                    avalue -= 360;
                }

                avalue = prompt("Enter angle (deg): ", avalue);

                if (avalue !== null) {
                    sket.con.type[  ncon] = "A";
                    sket.con.index1[ncon] = ibest;
                    sket.con.index2[ncon] = -1;
                    sket.con.value[ ncon] = avalue;

                    sket.pnt.lbl[ibest] += "A";
                }
            }

        // "h" --- Segment is horizontal
        } else if (myKeyPress == "h" || myKeyPress == "H") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sket.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sket.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (sket.seg.lbl[ibest].indexOf("I") >= 0) {
                alert("'I' constraint already on this Segment");
            } else if (Math.abs(sket.seg.dip[ibest]) > 1e-3) {
                alert("'H' constraint cannot be set on a Cirarc");
            } else {
                var ibeg  = sket.seg.ibeg[ibest];
                var iend  = sket.seg.iend[ibest];

                sket.con.type[  ncon] = "H";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[ibest] += "H";
            }

        // "v" --- Segment is vertical
        } else if (myKeyPress == "v" || myKeyPress == "V") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sket.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (sket.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sket.seg.lbl[ibest].indexOf("I") >= 0) {
                alert("'I' constraint already on this Segment");
            } else if (Math.abs(sket.seg.dip[ibest]) > 1e-3) {
                alert("'V' constraint cannot be set on a Cirarc");
            } else {
                var ibeg  = sket.seg.ibeg[ibest];
                var iend  = sket.seg.iend[ibest];

                sket.con.type[  ncon] = "V";
                sket.con.index1[ncon] = ibeg;
                sket.con.index2[ncon] = iend;
                sket.con.value[ ncon] = 0;

                sket.seg.lbl[ibest] += "V";
            }

        // "i" --- inclination of Segment
        } else if (myKeyPress == "i" || myKeyPress == "I") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            var ivalue;
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sket.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sket.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (Math.abs(sket.seg.dip[ibest]) > 1e-3) {
                alert("'I' constraint cannot be set on a Cirarc");
            } else if (sket.seg.lbl[ibest].indexOf("I") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[ icon] == "I" && sket.con.index1[icon] == ibest) {
                        ivalue = prompt("Enter inclination angle (deg): ", sket.con.value[icon]);
                        if (ivalue !== null) {
                            sket.con.value[icon] = ivalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg  = sket.seg.ibeg[ibest];
                var iend  = sket.seg.iend[ibest];

                var iprompt = Math.atan2(sket.pnt.y[ibeg]-sket.pnt.y[iend], sket.pnt.x[iend]-sket.pnt.x[ibeg]) * 180/Math.PI;

                ivalue = prompt("Enter inclination angle (deg): ", iprompt);
                if (ivalue !== null) {
                    sket.con.type[  ncon] = "I";
                    sket.con.index1[ncon] = ibeg;
                    sket.con.index2[ncon] = iend;
                    sket.con.value[ ncon] = ivalue;

                    sket.seg.lbl[ibest] += "I";
                }
            }

        // "l" --- length of Segment
        } else if (myKeyPress == "l" || myKeyPress == "L") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            var lvalue;
            if (ibest < 0) {
                alert("No Segment found");
            } else if (Math.abs(sket.seg.dip[ibest]) > 1e-3) {
                alert("'L' constraint cannot be set on a Cirarc");
            } else if (sket.seg.lbl[ibest].indexOf("L") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[ icon] == "L" && sket.con.index1[icon] == ibest) {
                        lvalue = prompt("Enter length: ", sket.con.value[icon]);
                        if (lvalue !== null) {
                            sket.con.value[icon] = lvalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg  = sket.seg.ibeg[ibest];
                var iend  = sket.seg.iend[ibest];

                var xbeg = sket.pnt.x[ibeg];
                var ybeg = sket.pnt.y[ibeg];
                var xend = sket.pnt.x[iend];
                var yend = sket.pnt.y[iend];
                var len  = Math.sqrt(Math.pow(xend-xbeg, 2) + Math.pow(yend-ybeg, 2));

                if (sket.scale === undefined) {
                    lvalue = prompt("Enter length: ");
                    if (lvalue !== null) {
                        sket.scale = lvalue / len;
                    }
                } else {
                    lvalue = prompt("Enter length: ", sket.scale*len);
                }

                if (lvalue !== null) {
                    sket.con.type[  ncon] = "L";
                    sket.con.index1[ncon] = ibeg;
                    sket.con.index2[ncon] = iend;
                    sket.con.value[ ncon] = lvalue;

                    sket.seg.lbl[ibest] += "L";
                }
            }

        // "r" --- radius of Segment
        } else if (myKeyPress == "r" || myKeyPress == "R") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            var rvalue;
            if (ibest < 0) {
                alert("No Cirarc found");
            } else if (Math.abs(sket.seg.dip[ibest]) < 1e-3) {
                alert("'R' constraint cannot be set for a Linseg");
            } else if (sket.seg.lbl[ibest].indexOf("R") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[ icon] == "R" && sket.con.index1[icon] == ibest) {
                        rvalue = prompt("Enter radius: ", sket.con.value[icon]);
                        if (rvalue !== null) {
                            sket.con.value[icon] = rvalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = sket.seg.ibeg[ibest];
                var iend = sket.seg.iend[ibest];

                var xa = sket.pnt.x[ibeg];
                var ya = sket.pnt.y[ibeg];
                var xb = sket.pnt.x[iend];
                var yb = sket.pnt.y[iend];

                var d  = sket.seg.dip[ibest];

                var L   = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));
                var rad = (L*L + 4*d*d) / (8*d);

                if (sket.scale === undefined) {
                    if (Number(sket.seg.dip[ibest]) >= 0) {
                        rvalue = prompt("Enter radius (should be positive as drawn): ");
                    } else {
                        rvalue = prompt("Enter radius (should be negative as drawn): ");
                    }
                    if (rvalue !== null) {
                        sket.scale = rvalue / rad;
                    }
                } else {
                    rvalue = prompt("Enter radius: ", sket.scale*rad);
                }

                if (rvalue !== null) {
                    sket.con.type[  ncon] = "R";
                    sket.con.index1[ncon] = ibeg;
                    sket.con.index2[ncon] = iend;
                    sket.con.value[ ncon] = rvalue;

                    sket.seg.lbl[ibest] += "R";
                }
            }

        // "s" --- sweep angle of Segment
        } else if (myKeyPress == "s" || myKeyPress == "S") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchSegment();
            var avalue;
            if (ibest < 0) {
                alert("No Cirarc found");
            } else if (Math.abs(sket.seg.dip[ibest]) < 1e-3) {
                alert("'S' constraint cannot be set for a Linseg");
            } else if (sket.seg.lbl[ibest].indexOf("S") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[icon] == "S" && sket.con.index1[icon] == ibest) {
                        avalue = prompt("Enter sweep angle (deg): ", sket.con.value[icon]);
                        if (avalue !== null) {
                            sket.con.value[icon] = avalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = sket.seg.ibeg[ibest];
                var iend = sket.seg.iend[ibest];

                if (Number(sket.seg.dip[ibest]) >= 0) {
                    avalue = prompt("Enter sweep angle (deg, should be positive as drawn): ");
                } else {
                    avalue = prompt("Enter sweep angle (deg, should be negative as drawn): ");
                }
                if (avalue !== null) {
                    sket.con.type[  ncon] = "S";
                    sket.con.index1[ncon] = ibeg;
                    sket.con.index2[ncon] = iend;
                    sket.con.value[ ncon] = avalue;

                    sket.seg.lbl[ibest] += "S";
                }
            }

        // "w" --- width (dx) between Points
        } else if (myKeyPress == "w" || myKeyPress == "W") {
            var bestPnt = getClosestSketchPoint();
            var bestCon = getClosestSketchConstraint();

            if (bestCon !== undefined) {
                var xmid = (  sket.pnt.x[sket.con.index1[bestCon]]
                            + sket.pnt.x[sket.con.index2[bestCon]]) / 2;
                var ymid = (2*sket.pnt.y[sket.con.index1[bestCon]]
                            + sket.pnt.y[sket.con.index2[bestCon]]) / 3;

                var dpnt = Math.abs(sket.pnt.x[bestPnt] - wv.cursorX)
                         + Math.abs(sket.pnt.y[bestPnt] - wv.cursorY);
                var dcon = Math.abs(xmid                - wv.cursorX)
                         + Math.abs(ymid                - wv.cursorY);

                if (dcon < dpnt) {
                    var wvalue = prompt("Enter width: ", sket.con.value[bestCon]);
                    if (wvalue !== null) {
                        sket.con.value[bestCon] = wvalue;
                    }

                    sket.bestPoint = undefined;
                    sket.mode      = 3;
                    sket.suggest   = undefined;
                } else {
                    sket.basePoint = bestPnt;
                    sket.mode      = 4;
                    sket.suggest   = undefined;
                }
            } else {
                sket.basePoint = bestPnt;
                sket.mode      = 4;
                sket.suggest   = undefined;
            }

        // "d" --- depth (dy) between Points
        } else if (myKeyPress == "d" || myKeyPress == "D") {
            var bestPnt = getClosestSketchPoint();
            var bestCon = getClosestSketchConstraint();

            if (bestCon !== undefined) {
                var xmid = (2*sket.pnt.x[sket.con.index1[bestCon]]
                            + sket.pnt.x[sket.con.index2[bestCon]]) / 3;
                var ymid = (  sket.pnt.y[sket.con.index1[bestCon]]
                            + sket.pnt.y[sket.con.index2[bestCon]]) / 2;

                var dpnt = Math.abs(sket.pnt.x[bestPnt] - wv.cursorX)
                         + Math.abs(sket.pnt.y[bestPnt] - wv.cursorY);
                var dcon = Math.abs(xmid                - wv.cursorX)
                         + Math.abs(ymid                - wv.cursorY);

                if (dcon < dpnt) {
                    var dvalue = prompt("Enter depth: ", sket.con.value[bestCon]);
                    if (dvalue !== null) {
                        sket.con.value[bestCon] = dvalue;
                    }

                    sket.bestPoint = undefined;
                    sket.mode      = 3;
                    sket.suggest   = undefined;
                } else {
                    sket.basePoint = bestPnt;
                    sket.mode      = 5;
                    sket.suggest   = undefined;
                }
            } else {
                sket.basePoint = bestPnt;
                sket.mode      = 5;
                sket.suggest   = undefined;
            }

        // "<" --- remove constraint(s) from nearset Point or Segment
        } else if (myKeyPress == "<") {
            sket.mode    = 3;
            sket.suggest = undefined;

            var ibest = getClosestSketchPoint();
            var jbest = getClosestSketchSegment();
            var kbest = getClosestSketchConstraint();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sket.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sket.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sket.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sket.seg.ym[jbest]);
            }
            var dkbest = 999999;
            if (kbest !== undefined) {
                if        (sket.con.type[kbest] == "W") {
                    var xmid = (  sket.pnt.x[sket.con.index1[kbest]]
                                + sket.pnt.x[sket.con.index2[kbest]]) / 2;
                    var ymid = (2*sket.pnt.y[sket.con.index1[kbest]]
                                + sket.pnt.y[sket.con.index2[kbest]]) / 3;

                    dkbest = Math.abs(wv.cursorX - xmid)
                           + Math.abs(wv.cursorY - ymid);
                } else if (sket.con.type[kbest] == "D") {
                    var xmid = (2*sket.pnt.x[sket.con.index1[kbest]]
                                + sket.pnt.x[sket.con.index2[kbest]]) / 3;
                    var ymid = (  sket.pnt.y[sket.con.index1[kbest]]
                                + sket.pnt.y[sket.con.index2[kbest]]) / 2;

                    dkbest = Math.abs(wv.cursorX - xmid)
                           + Math.abs(wv.cursorY - ymid);
                }
            }

            // remove (selected) constraint at Point ibest
            if (ibest >= 0 && dibest < djbest && dibest < dkbest) {
                if (sket.pnt.lbl[ibest].length == 0) {
                    postMessage("no constraint at Point "+ibest);
                } else if (sket.pnt.lbl[ibest] == "Z") {
                    postMessage("cannot delete a 'Z' constraint");
                } else {
                    var jtype;
                    if (sket.pnt.lbl[ibest].length == 1) {
                        jtype = sket.pnt.lbl[ibest];
                    } else {
                        jtype = prompt("More than one constraint at Point "+ibest+
                                      "\nEnter one character from \""+sket.pnt.lbl[ibest]+
                                      "\" to remove");
                        if (jtype === null) {
                            jtype = "&";
                        }
                    }
                    jtype = jtype.toUpperCase();

                    if        (ibest == 0 && sket.relative == 1 && jtype == "X") {
                        alert("X cannot be deleted from first point if in relative mode");
                    } else if (ibest == 0 && sket.relative == 1 && jtype == "Y") {
                        alert("Y cannot be deleted from first point if in relative mode");
                    } else {
                        for (var icon = 0; icon < sket.con.type.length; icon++) {
                            if (sket.con.index1[icon] == ibest &&
                                sket.con.index2[icon] <  0     &&
                                sket.con.type[  icon] == jtype   ) {
                                postMessage("removing \""+jtype+"\" contraint from Point "+ibest);

                                sket.con.type.splice(  icon, 1);
                                sket.con.index1.splice(icon, 1);
                                sket.con.index2.splice(icon, 1);
                                sket.con.value.splice( icon, 1);

                                sket.pnt.lbl[ibest] = sket.pnt.lbl[ibest].split(jtype).join("");
                                break;
                            }
                        }
                    }
                }

            // remove constraints at Segment jbest
            } else if (jbest >= 0 && djbest < dkbest) {
                if (sket.seg.lbl[jbest].length == 0) {
                    postMessage("no constraint at Segment "+jbest);
                } else if (sket.seg.lbl[ibest] == "Z") {
                    postMessage("cannot delete a 'Z' constraint");
                } else {
                    var jtype;
                    if (sket.seg.lbl[jbest].length == 1) {
                        jtype = sket.seg.lbl[jbest];
                    } else {
                        jtype = prompt("More than one constraint at Segment "+jbest+
                                      "\nEnter one character from \""+sket.seg.lbl[jbest]+
                                      "\" to remove");
                        if (jtype === null) {
                            jtype = "&";
                        }
                    }
                    jtype = jtype.toUpperCase();

                    for (var icon = 0; icon < sket.con.type.length; icon++) {
                        if (sket.con.index1[icon] == jbest &&
                            sket.con.index2[icon] >= 0     &&
                            sket.con.type[  icon] == jtype   ) {
                            postMessage("removing \""+jtype+"\" constraint from Segment "+jbest);

                            sket.con.type.splice(  icon, 1);
                            sket.con.index1.splice(icon, 1);
                            sket.con.index2.splice(icon, 1);
                            sket.con.value.splice( icon, 1);

                            sket.seg.lbl[jbest] = sket.seg.lbl[jbest].split(jtype).join("");
                            break;
                        }
                    }
                }

            // remove H or W constraint
            } else {
                postMessage("removing constraint "+kbest);

                sket.con.type.splice(  kbest, 1);
                sket.con.index1.splice(kbest, 1);
                sket.con.index2.splice(kbest, 1);
                sket.con.value.splice( kbest, 1);
            }

        // "?" --- examine Point or Segment
        } else if (myKeyPress == "?") {
            var ibest = getClosestSketchPoint();
            var jbest = getClosestSketchSegment();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sket.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sket.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sket.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sket.seg.ym[jbest]);
            }

            // examine Point ibest
            if (ibest >= 0 && dibest < djbest) {
                if (sket.pnt.lbl[ibest].length == 0) {
                    postMessage("Point "+(ibest+1)+" has no constraints");
                } else {
                    postMessage("Point "+(ibest+1)+" has constraints: "+sket.pnt.lbl[ibest]);
                }

            // examine Segment jbest
            } else if (jbest >= 0) {
                if        (sket.seg.type[jbest] == "L") {
                    var segtype = "linseg";
                } else if (sket.seg.type[jbest] == "C") {
                    var segtype = "cirarc";
                } else if (sket.seg.type[jbest] == "S") {
                    var segtype = "spline";
                } else if (sket.seg.type[jbest] == "B") {
                    var segtype = "bezier";
                } else {
                    var segtype = "**unknown type**";
                }
                if (sket.seg.lbl[jbest].length == 0) {
                    postMessage("Segment "+(jbest+1)+" ("+segtype+") has no constraints");
                } else {
                    postMessage("Segment "+(jbest+1)+" ("+segtype+") has constraints: "+sket.seg.lbl[jbest]);
                }
            }

        // "@" --- dump sket structure
        } else if (myKeyPress == "@") {
            postMessage("sket.pnt->"+JSON.stringify(sket.pnt));
            postMessage("sket.seg->"+JSON.stringify(sket.seg));
            postMessage("sket.var->"+JSON.stringify(sket.var));
            postMessage("sket.con->"+JSON.stringify(sket.con));

        // invalid option
        } else {
            var ibest = getClosestSketchPoint();
            var jbest = getClosestSketchSegment();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sket.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sket.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sket.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sket.seg.ym[jbest]);
            }

            // remove constraints at Point ibest
            if (ibest >= 0 && dibest < djbest) {
                postMessage("Valid options are 'x', 'y', 'p', 't', 'w', 'd', and '<'");
            } else {
                postMessage("Valid options are 'h', 'v', 'i', 'l', 'r', 's', and '<'");
            }
        }

        // if .scale has been set, try to update .xorig or .yorig
        //    if a suitable constraint exists
        if (sket.scale !== undefined) {
            if (sket.xorig === undefined) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[icon] == "X") {
                        sket.xorig = sket.pnt.x[sket.con.index1[icon]] - sket.con.value[icon] / sket.scale;
                        break;
                    }
                }
            }
            if (sket.yorig === undefined) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sket.con.type[icon] == "Y") {
                        sket.yorig = sket.pnt.y[sket.con.index1[icon]] + sket.con.value[icon] / sket.scale;
                        break;
                    }
                }
            }
        }

    // sket.mode 4 options: w
    } else if (sket.mode == 4) {
        sket.mode    = 3;
        sket.suggest = undefined;

        var target = getClosestSketchPoint();

        if (target == sket.basePoint) {
            alert("Width not set since base and target Point are the same");
        } else {
            var wvalue;
            if (Number(sket.pnt.x[sket.basePoint]) <= Number(sket.pnt.x[target])) {
                wvalue = prompt("Enter width (should be positive as drawn): ", wvalue);
            } else {
                wvalue = prompt("Enter width (should be negative as drawn): ", wvalue);
            }
            if (wvalue !== null) {
                if (isNaN(wvalue) !== true) {
                    if (sket.pnt.x[target] > sket.pnt.x[sket.basePoint]) {
                        wvalue = +Math.abs(wvalue);
                    } else {
                        wvalue = -Math.abs(wvalue);
                    }
                }

                sket.con.type[  ncon] = "W";
                sket.con.index1[ncon] = sket.basePoint;
                sket.con.index2[ncon] = target;
                sket.con.value[ ncon] = wvalue;
            }
        }

        sket.basePoint = undefined;

    // sket.mode 5 options: d
    } else if (sket.mode == 5) {
        sket.mode    = 3;
        sket.suggest = undefined;

        var target = getClosestSketchPoint();

        if (target == sket.basePoint) {
            alert("Depth not set since base and target Point are the same");
        } else {
            var dvalue;
            if (Number(sket.pnt.y[sket.basePoint]) >= Number(sket.pnt.y[target])) {
                dvalue = prompt("Enter depth (should be positive as drawn): ", dvalue);
            } else {
                dvalue = prompt("Enter depth (should be negative as drawn): ", dvalue);
            }
            if (dvalue !== null) {
                if (isNaN(dvalue) !== true) {
                    if (sket.pnt.y[target] > sket.pnt.y[sket.basePoint]) {
                        dvalue = +Math.abs(dvalue);
                    } else {
                        dvalue = -Math.abs(dvalue);
                    }
                }

                sket.con.type[  ncon] = "D";
                sket.con.index1[ncon] = sket.basePoint;
                sket.con.index2[ncon] = target;
                sket.con.value[ ncon] = dvalue;
            }
        }

        sket.basePoint = undefined;
    }

    drawSketch();
}


//
// called when "solveButton" is pressed (defined in ESP-sketch.html)
//
function solveSketchPre()
{
    // alert("in solveSketchPre()");

    var npnt = sket.pnt.x.length;
    var nvar = sket.var.name.length;
    var ncon = sket.con.type.length;
    var nseg = sket.seg.type.length;

    // if (nvar > ncon) {
    //     alert("Sketch is under-constrained since nvar="+nvar+" but ncon="+ncon);
    //     return;
    // } else if (nvar < ncon) {
    //     alert("Sketch is over-constrained since nvar="+nvar+" but ncon="+ncon);
    //     return;
    // }

    if (sket.xorig === undefined || isNaN(sket.xorig) ||
        sket.yorig === undefined || isNaN(sket.yorig) ||
        sket.scale === undefined || isNaN(sket.scale)   ) {

        // get the size of the canvas
        var canvas = document.getElementById("sketcher");
        var width  = canvas.clientWidth;
        var height = canvas.clientHeight;

        // guess that scale of body is "5" and that it is centered at origin
        sket.scale = 5 / width;
        sket.xorig = width  / 2;
        sket.yorig = height / 2;
    }

    // create a message to send to the server
    var message = "solveSketch|";

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        if (ipnt == 0) {
            var im1 = npnt - 1;
        } else {
            var im1 = ipnt - 1;
        }

        message += Number(sket.scale * (sket.pnt.x[ipnt] - sket.xorig)).toFixed(6) + ";";
        message += Number(sket.scale * (sket.yorig - sket.pnt.y[ipnt])).toFixed(6) + ";";

        if (ipnt > 0) {
            message += Number(sket.scale *  sket.seg.dip[ipnt-1]      ).toFixed(6) + ";";
        } else if (sket.seg.iend[nseg-1] == 0) {
            message += Number(sket.scale *  sket.seg.dip[npnt-1]      ).toFixed(6) + ";";
        } else {
            message += "0;";
        }
    }

    message += "|";

    for (var icon = 0; icon < ncon; icon++) {
        message += sket.con.type[icon] + ";" + (sket.con.index1[icon]+1) + ";";
        if (sket.con.index2[icon] >= 0) {
            message += (sket.con.index2[icon] + 1) + ";" + sket.con.value[icon] + ";";
        } else {
            message +=  sket.con.index2[icon]      + ";" + sket.con.value[icon] + ";";
        }
    }

    message += "|";

    browserToServer(message);
}


//
// update the variables associated with a Sketch
//
function solveSketchPost(text)
{
    // alert("solveSketchPost(text="+text+")");

    // break up based upon the "|"
    var textList = text.split("|");

    if        (textList[1].substring(0,7) == "ERROR::") {
        postMessage(textList[1]);
        alert("Sketch was not successfully solved, but no hints are available\n" +
              "Note: try adding a R or S constraint to a circular arc.");
        return;
    } else if (textList[1].substring(0,4) == "*del") {
        sket.suggest = textList[1];
        drawSketch();
        postMessage("Delete one of the constraints in red (using < key)\n" +
                    "Note: others may be possible as well.");
        return;
    } else if (textList[1].substring(0,4) == "*add") {
        sket.suggest = textList[1];
        drawSketch();
        postMessage("Add one of the constraints in green\n" +
                    "Note: others may be possible as well.");
        return;
    }

    var npnt = 0;
    for (var i = 0; i < textList[1].length; i++) {
        if (textList[1].charAt(i) == ';') {
            npnt++;
        }
    }
    npnt /= 3;

    // update the values (which are separated by semicolons)
    var values = textList[1].split(";");

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        if (ipnt == 0) {
            var im1 = npnt - 1;
        } else {
            var im1 = ipnt - 1;
        }

        sket.pnt.x[  ipnt] = Math.floor(sket.xorig + Number(values[3*ipnt  ]) / sket.scale);
        sket.pnt.y[  ipnt] = Math.floor(sket.yorig - Number(values[3*ipnt+1]) / sket.scale);
        sket.seg.dip[im1 ] = Math.floor(             Number(values[3*ipnt+2]) / sket.scale);
    }

    // update data associated with segments
    for (var iseg = 0; iseg < sket.seg.type.length; iseg++) {
        updateSegData(iseg);
    }

    // set the mode
    sket.mode = 6;

    // redraw the Sketch with the updated values
    drawSketch();
}


//
// called when "quit" is pressed (defined in ESP-sketch.html)
//
function quitSketch()
{
    // alert("in quitSketch()");

    changeMode(0);
}


//
// called when "save" is pressed (defined in ESP-sketch.html)
//
function saveSketch()
{
    // alert("in saveSketch()");

    // if the scaling is not set, make sure that at least one length is given
    //    (or else NaN will be sent back to .csm)
    if        (sket.scale === undefined || isNaN(sket.scale)) {
        alert("At lest one length or radius must be set before saving");
        return;
    } else if (sket.xorig === undefined || isNaN(sket.xorig)) {
        alert("At least one \"X\" must be set before saving");
        return;
    } else if (sket.yorig === undefined || isNaN(sket.yorig)) {
        alert("At least one \"Y\" must be set before saving");
        return;
    }

    // alert user that an unsolved Sketch may cause .csm to not solve
    if (sket.mode != 6) {
        if (confirm("Sketch is not solved and may break build.  Continue?") !== true) {
            return;
        }
    }

    var npnt = sket.pnt.x.length;
    var nseg = sket.seg.type.length;
    var ncon = sket.con.type.length;

    // build a string of the variables
    var vars = "";
    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        vars += Number(sket.scale * (sket.pnt.x[ipnt] - sket.xorig)).toFixed(6) + ";";
        vars += Number(sket.scale * (sket.yorig - sket.pnt.y[ipnt])).toFixed(6) + ";";

        if (ipnt > 0) {
            vars += Number(sket.scale *  sket.seg.dip[ipnt-1]      ).toFixed(6) + ";";
        } else if (sket.seg.iend[nseg-1] == 0) {
            vars += Number(sket.scale *  sket.seg.dip[npnt-1]      ).toFixed(6) + ";";
        } else {
            vars += "0;";
        }
    }

    // build a string of the constraints
    var cons = "";
    for (var icon = 0; icon < ncon; icon++) {
        cons += sket.con.type[icon] + ";" + (sket.con.index1[icon]+1) + ";";
        if (sket.con.index2[icon] >= 0) {
            cons += (sket.con.index2[icon] + 1) + ";" + sket.con.value[icon] + ";";
        } else {
            cons +=  sket.con.index2[icon]      + ";" + sket.con.value[icon] + ";";
        }
    }
    cons = cons.replace(/\s/g, "");

    // build a string of the segments
    var segs = "";
    for (var iseg = 0; iseg < nseg; iseg++) {
        segs += sket.seg.type[iseg] + ";" + (sket.seg.ibeg[iseg]+1) + ";" +
                                            (sket.seg.iend[iseg]+1) + ";";
    }

    // create the message to be sent to the server
    var mesg = ""+sket.ibrch+"|"+vars+"|"+cons+"|"+segs+"|";

    // because of an apparent limit on the size of text
    //    messages that can be sent from the browser to the
    //    server, we need to send the new file back in
    //    pieces and then reassemble on the server
    var maxMessageSize = 1000;

    var ichar = 0;
    var part  = mesg.substring(ichar, ichar+maxMessageSize);
    browserToServer("saveSketchBeg|"+part);
    ichar += maxMessageSize;

    while (ichar < mesg.length) {
        part = mesg.substring(ichar, ichar+maxMessageSize);
        browserToServer("saveSketchMid|"+part);
        ichar += maxMessageSize;
    }

    browserToServer("saveSketchEnd|");
}


//
// undo last Sketch command
//
function undoSketch()
{
    // alert("undoSketch()");

    var iundo = (sket.nundo-1) % sket.maxundo;

    // if undo information does not exist, inform the user
    if (sket.undo[iundo]     === undefined ||
        sket.undo[iundo].pnt === undefined ||
        sket.undo[iundo].seg === undefined ||
        sket.undo[iundo].var === undefined ||
        sket.undo[iundo].con === undefined   ) {
        alert("There is nothing to undo");
        return;
    }

    // remove old Sketch info and initialize the arrays
    sket.pnt        = {};
    sket.pnt.x      = [];
    sket.pnt.y      = [];
    sket.pnt.lbl    = [];

    sket.seg        = {};
    sket.seg.type   = [];
    sket.seg.ibeg   = [];
    sket.seg.iend   = [];
    sket.seg.xm     = [];
    sket.seg.ym     = [];
    sket.seg.lbl    = [];

    sket.seg.dip    = [];
    sket.seg.xc     = [];
    sket.seg.yc     = [];
    sket.seg.rad    = [];
    sket.seg.tbeg   = [];
    sket.seg.tend   = [];

    sket.var        = {};
    sket.var.name   = [];
    sket.var.value  = [];

    sket.con        = {};
    sket.con.type   = [];
    sket.con.index1 = [];
    sket.con.index2 = [];
    sket.con.value  = [];

    // copy the undo information into Sketch information
    sket.mode = sket.undo[iundo].mode;

    for (var ipnt = 0; ipnt < sket.undo[iundo].pnt.x.length; ipnt++) {
        sket.pnt.x[  ipnt] =  sket.undo[iundo].pnt.x[  ipnt];
        sket.pnt.y[  ipnt] =  sket.undo[iundo].pnt.y[  ipnt];
        sket.pnt.lbl[ipnt] =  sket.undo[iundo].pnt.lbl[ipnt];
    }
    for (var iseg = 0; iseg < sket.undo[iundo].seg.type.length; iseg++) {
        sket.seg.type[iseg] = sket.undo[iundo].seg.type[iseg];
        sket.seg.ibeg[iseg] = sket.undo[iundo].seg.ibeg[iseg];
        sket.seg.iend[iseg] = sket.undo[iundo].seg.iend[iseg];
        sket.seg.xm[  iseg] = sket.undo[iundo].seg.xm[  iseg];
        sket.seg.ym[  iseg] = sket.undo[iundo].seg.ym[  iseg];
        sket.seg.lbl[ iseg] = sket.undo[iundo].seg.lbl[ iseg];
        sket.seg.dip[ iseg] = sket.undo[iundo].seg.dip[ iseg];
        sket.seg.xc[  iseg] = sket.undo[iundo].seg.xc[  iseg];
        sket.seg.yc[  iseg] = sket.undo[iundo].seg.yc[  iseg];
        sket.seg.rad[ iseg] = sket.undo[iundo].seg.rad[ iseg];
        sket.seg.tbeg[iseg] = sket.undo[iundo].seg.tbeg[iseg];
        sket.seg.tend[iseg] = sket.undo[iundo].seg.tend[iseg];
    }
    for (var ivar = 0; ivar <  sket.undo[iundo].var.name.length; ivar++) {
        sket.var.name[ ivar] = sket.undo[iundo].var.name[ ivar];
        sket.var.value[ivar] = sket.undo[iundo].var.value[ivar];
    }
    for (var icon = 0; icon <   sket.undo[iundo].con.type.length; icon++){
        sket.con.type[  icon] = sket.undo[iundo].con.type[  icon];
        sket.con.index1[icon] = sket.undo[iundo].con.index1[icon];
        sket.con.index2[icon] = sket.undo[iundo].con.index2[icon];
        sket.con.value[ icon] = sket.undo[iundo].con.value[ icon];
    }

    // remove the undo information
    sket.undo[iundo] = {};

    // decrement number of undos
    sket.nundo--;

    // update the display
    drawSketch();
}


//
// save an undo snapshot
//
function saveSketchUndo()
{
    // alert("saveSketchUndo()");

    var iundo = sket.nundo % sket.maxundo;

    // remove old undo info and initialize the arrays
    sket.undo[iundo]            = {};

    sket.undo[iundo].pnt        = {};
    sket.undo[iundo].pnt.x      = [];
    sket.undo[iundo].pnt.y      = [];
    sket.undo[iundo].pnt.lbl    = [];

    sket.undo[iundo].seg        = {};
    sket.undo[iundo].seg.type   = [];
    sket.undo[iundo].seg.ibeg   = [];
    sket.undo[iundo].seg.iend   = [];
    sket.undo[iundo].seg.xm     = [];
    sket.undo[iundo].seg.ym     = [];
    sket.undo[iundo].seg.lbl    = [];

    sket.undo[iundo].seg.dip    = [];
    sket.undo[iundo].seg.xc     = [];
    sket.undo[iundo].seg.yc     = [];
    sket.undo[iundo].seg.rad    = [];
    sket.undo[iundo].seg.tbeg   = [];
    sket.undo[iundo].seg.tend   = [];

    sket.undo[iundo].var        = {};
    sket.undo[iundo].var.name   = [];
    sket.undo[iundo].var.value  = [];

    sket.undo[iundo].con        = {};
    sket.undo[iundo].con.type   = [];
    sket.undo[iundo].con.index1 = [];
    sket.undo[iundo].con.index2 = [];
    sket.undo[iundo].con.value  = [];

    // copy the Sketch information into the undo storage
    sket.undo[iundo].mode = sket.mode;

    for (var ipnt = 0; ipnt < sket.pnt.x.length; ipnt++) {
        sket.undo[iundo].pnt.x[  ipnt] = sket.pnt.x[  ipnt];
        sket.undo[iundo].pnt.y[  ipnt] = sket.pnt.y[  ipnt];
        sket.undo[iundo].pnt.lbl[ipnt] = sket.pnt.lbl[ipnt];
    }
    for (var iseg = 0; iseg < sket.seg.type.length; iseg++) {
        sket.undo[iundo].seg.type[iseg] = sket.seg.type[iseg];
        sket.undo[iundo].seg.ibeg[iseg] = sket.seg.ibeg[iseg];
        sket.undo[iundo].seg.iend[iseg] = sket.seg.iend[iseg];
        sket.undo[iundo].seg.xm[  iseg] = sket.seg.xm[  iseg];
        sket.undo[iundo].seg.ym[  iseg] = sket.seg.ym[  iseg];
        sket.undo[iundo].seg.lbl[ iseg] = sket.seg.lbl[ iseg];
        sket.undo[iundo].seg.dip[ iseg] = sket.seg.dip[ iseg];
        sket.undo[iundo].seg.xc[  iseg] = sket.seg.xc[  iseg];
        sket.undo[iundo].seg.yc[  iseg] = sket.seg.yc[  iseg];
        sket.undo[iundo].seg.rad[ iseg] = sket.seg.rad[ iseg];
        sket.undo[iundo].seg.tbeg[iseg] = sket.seg.tbeg[iseg];
        sket.undo[iundo].seg.tend[iseg] = sket.seg.tend[iseg];
    }
    for (var ivar = 0; ivar < sket.var.name.length; ivar++) {
        sket.undo[iundo].var.name[ ivar] = sket.var.name[ ivar];
        sket.undo[iundo].var.value[ivar] = sket.var.value[ivar];
    }
    for (var icon = 0; icon < sket.con.type.length; icon++){
        sket.undo[iundo].con.type[  icon] = sket.con.type[  icon];
        sket.undo[iundo].con.index1[icon] = sket.con.index1[icon];
        sket.undo[iundo].con.index2[icon] = sket.con.index2[icon];
        sket.undo[iundo].con.value[ icon] = sket.con.value[ icon];
    }

    // increment number of undos
    sket.nundo++;
}


//
// update the variables associated with a segment
//
function updateSegData(iseg)
{
    // alert("updateSegData(iseg="+iseg+")");

    sket.seg.xm[iseg] = (sket.pnt.x[sket.seg.ibeg[iseg]] + sket.pnt.x[sket.seg.iend[iseg]]) / 2;
    sket.seg.ym[iseg] = (sket.pnt.y[sket.seg.ibeg[iseg]] + sket.pnt.y[sket.seg.iend[iseg]]) / 2;

    if (Math.abs(sket.seg.dip[iseg]) > 1.0e-6) {
        var xa = sket.pnt.x[sket.seg.ibeg[iseg]];
        var ya = sket.pnt.y[sket.seg.ibeg[iseg]];
        var xb = sket.pnt.x[sket.seg.iend[iseg]];
        var yb = sket.pnt.y[sket.seg.iend[iseg]];
        var d  =            sket.seg.dip[ iseg];

        var L  = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));
        var R  = (L*L + 4*d*d) / (8*d);

        sket.seg.xc[iseg] = (xa + xb) / 2 + (R - d) * (yb - ya) / L;
        sket.seg.yc[iseg] = (ya + yb) / 2 - (R - d) * (xb - xa) / L;

        sket.seg.rad[ iseg] = R;
        sket.seg.tbeg[iseg] = Math.atan2(ya-sket.seg.yc[iseg], xa-sket.seg.xc[iseg]);
        sket.seg.tend[iseg] = Math.atan2(yb-sket.seg.yc[iseg], xb-sket.seg.xc[iseg]);

        if (R > 0) {
            sket.seg.tbeg[iseg] -= 2 * Math.PI;
            while (sket.seg.tbeg[iseg] < sket.seg.tend[iseg]) {
                sket.seg.tbeg[iseg] += 2 * Math.PI;
            }
        } else {
            sket.seg.tbeg[iseg] += 2 * Math.PI;
            while (sket.seg.tbeg[iseg] > sket.seg.tend[iseg]) {
                sket.seg.tbeg[iseg] -= 2 * Math.PI;
            }
        }

        var tavg = (sket.seg.tbeg[iseg] + sket.seg.tend[iseg]) / 2;
        sket.seg.xm[iseg] = sket.seg.xc[iseg] + Math.cos(tavg) * Math.abs(R);
        sket.seg.ym[iseg] = sket.seg.yc[iseg] + Math.sin(tavg) * Math.abs(R);
    }
}


//
// Draw the current Sketch
//
function drawSketch()
{
    // alert("drawSketch()");

    // set up context
    var canvas  = document.getElementById("sketcher");
    var context = canvas.getContext("2d");

    // clear the screen
    context.clearRect(0, 0, canvas.width, canvas.height);

    context.lineCap  = "round";
    context.lineJoin = "round";

    // determine number of Points and Segments in Sketch
    var npnt = sket.pnt.x.length;
    var nseg = sket.seg.type.length;
    var ncon = sket.con.type.length;

    // return if no Points
    if (npnt <= 0) {
        var button = document.getElementById("buildButton");
        button.style.backgroundColor =  "#FFFF3F";
        button["innerHTML"] = "Initializing...";

        return;
    }

    // axes
    context.strokeStyle = "yellow";
    context.lineWidth   = 1;

    context.beginPath();
    context.moveTo(sket.pnt.x[0], -100);
    context.lineTo(sket.pnt.x[0], 1000);
    context.stroke();

    context.beginPath();
    context.moveTo(-100, sket.pnt.y[0]);
    context.lineTo(1000, sket.pnt.y[0]);
    context.stroke();

    // draw the part of the Sketch that is complete
    if (nseg > 0) {

        // we are editting the Cirarc at the end, so it has to
        //    be drawn separately (since it will be drawn in red)
        if (sket.mode == 2) {

            // the part of the Sketch that is done
            context.strokeStyle = "black";
            context.lineWidth   = 3;
            context.beginPath();

            var ibeg = sket.seg.ibeg[0];
            var iend = sket.seg.iend[0];
            context.moveTo(sket.pnt.x[ibeg], sket.pnt.y[ibeg]);

            for (var iseg = 0; iseg < nseg-1; iseg++) {
                iend = sket.seg.iend[iseg];

                if (Math.abs(sket.seg.dip[iseg]) < 1e-6) {
                    context.lineTo(sket.pnt.x[iend], sket.pnt.y[iend]);
                } else if (sket.seg.rad[iseg] > 0) {
                    context.arc(sket.seg.xc[iseg], sket.seg.yc[iseg],
                                +sket.seg.rad[iseg], sket.seg.tbeg[iseg], sket.seg.tend[iseg], 1);
                } else {
                    context.arc(sket.seg.xc[iseg], sket.seg.yc[iseg],
                                -sket.seg.rad[iseg], sket.seg.tbeg[iseg], sket.seg.tend[iseg], 0);
                }
            }

            context.stroke();

            // the Cirarc at the end which is to be rendered in red
            context.strokeStyle = "red";
            context.beginPath();

            if (sket.seg.dip[nseg-1] > 0) {
                context.arc(sket.seg.xc[nseg-1], sket.seg.yc[nseg-1],
                            +sket.seg.rad[nseg-1], sket.seg.tbeg[nseg-1], sket.seg.tend[nseg-1], 1);
            } else {
                context.arc(sket.seg.xc[nseg-1], sket.seg.yc[nseg-1],
                            -sket.seg.rad[nseg-1], sket.seg.tbeg[nseg-1], sket.seg.tend[nseg-1], 0);
            }

            context.stroke();

        // we are NOT editting the Cirarc at the end
        } else {
            context.strokeStyle = "black";
            context.lineWidth   = 3;
            context.beginPath();

            var ibeg = sket.seg.ibeg[0];
            var iend = sket.seg.iend[0];
            context.moveTo(sket.pnt.x[ibeg], sket.pnt.y[ibeg]);

            for (var iseg = 0; iseg < nseg; iseg++) {
                iend = sket.seg.iend[iseg];

                if (sket.seg.type[iseg] == "B") {
                    var Xbezier = [];
                    var Ybezier = [];
                    var nbezier = 2;

                    Xbezier[0] = sket.pnt.x[ibeg];
                    Ybezier[0] = sket.pnt.y[ibeg];
                    Xbezier[1] = sket.pnt.x[iend];
                    Ybezier[1] = sket.pnt.y[iend];

                    while (iseg < nseg-1) {
                        if (sket.seg.type[iseg+1] == "B") {
                            iseg++;
                            iend = sket.seg.iend[iseg];
                            Xbezier[nbezier] = sket.pnt.x[iend];
                            Ybezier[nbezier] = sket.pnt.y[iend];
                            nbezier++;
                        } else {
                            break;
                        }
                    }

                    drawBezier(context, Xbezier, Ybezier);
//$$$             } else if (sket.seg.type[iseg] == "S") {
//$$$                 var Xspline = [];
//$$$                 var Yspline = [];
//$$$                 var nspline = 2;
//
//$$$                 Xspline[0] = sket.pnt.x[ibeg];
//$$$                 Yspline[0] = sket.pnt.y[ibeg];
//$$$                 Xspline[1] = sket.pnt.x[iend];
//$$$                 Yspline[1] = sket.pnt.y[iend];
//
//$$$                 while (iseg < nseg-1) {
//$$$                     if (sket.seg.type[iseg+1] == "S") {
//$$$                         iseg++;
//$$$                         iend = sket.seg.iend[iseg];
//$$$                         Xspline[nspline] = sket.pnt.x[iend];
//$$$                         Yspline[nspline] = sket.pnt.y[iend];
//$$$                         nspline++;
//$$$                     } else {
//$$$                         break;
//$$$                     }
//$$$                 }
//
//$$$                 drawSpline(context, Xspline, Yspline);
                } else if (Math.abs(sket.seg.dip[iseg]) < 1e-6) {
                    context.lineTo(sket.pnt.x[iend], sket.pnt.y[iend]);
                } else if (sket.seg.dip[iseg] > 0) {
                    context.arc(sket.seg.xc[iseg], sket.seg.yc[iseg],
                                +sket.seg.rad[iseg], sket.seg.tbeg[iseg], sket.seg.tend[iseg], 1);
                } else {
                    context.arc(sket.seg.xc[iseg], sket.seg.yc[iseg],
                                -sket.seg.rad[iseg], sket.seg.tbeg[iseg], sket.seg.tend[iseg], 0);
                }

                ibeg = iend;
            }

            context.stroke();

            // if Sketch is closed, fill it
            if (iend == sket.seg.ibeg[0]) {
                if (sket.var.name.length == sket.con.type.length) {
                    context.fillStyle = "#c0f0c0";  // light green
                    context.fill();
                } else {
                    context.fillStyle = "#e0e0e0";  // grey
                    context.fill();
                }
            }

            // if draw straight lines for Segments associated with splines
            //    and beziers
            context.strokeStyle = "cyan";
            context.lineWidth   = 1;
            for (iseg = 0; iseg < nseg; iseg++) {
                if (sket.seg.type[iseg] == "B" ||
                    sket.seg.type[iseg] == "S"   ) {

                    ibeg = sket.seg.ibeg[iseg];
                    iend = sket.seg.iend[iseg];
                    context.beginPath();
                    context.moveTo(sket.pnt.x[ibeg], sket.pnt.y[ibeg]);
                    context.lineTo(sket.pnt.x[iend], sket.pnt.y[iend]);
                    context.stroke();
                }
            }

            // we are creating W constraint
            if (sket.mode == 4) {
                var xmid = (  sket.pnt.x[sket.basePoint] + wv.cursorX) / 2;
                var ymid = (2*sket.pnt.y[sket.basePoint] + wv.cursorY) / 3;

                context.strokeStyle = "red";
                context.lineWidth   = 1;
                context.beginPath();
                context.moveTo(sket.pnt.x[sket.basePoint], sket.pnt.y[sket.basePoint]);
                context.lineTo(sket.pnt.x[sket.basePoint], ymid                      );
                context.lineTo(wv.cursorX,                 ymid                      );
                context.lineTo(wv.cursorX,                 wv.cursorY                );
                context.stroke();
            }

            // we are creating D constraint
            if (sket.mode == 5) {
                var xmid = (2*sket.pnt.x[sket.basePoint] + wv.cursorX) / 3;
                var ymid = (  sket.pnt.y[sket.basePoint] + wv.cursorY) / 2;

                context.strokeStyle = "red";
                context.lineWidth   = 1;
                context.beginPath();
                context.moveTo(sket.pnt.x[sket.basePoint], sket.pnt.y[sket.basePoint]);
                context.lineTo(xmid,                       sket.pnt.y[sket.basePoint]);
                context.lineTo(xmid,                       wv.cursorY                );
                context.lineTo(wv.cursorX,                 wv.cursorY                );
                context.stroke();
            }
        }
    }

    // draw the Points
    context.fillStyle = "black";
    context.fillRect(sket.pnt.x[0]-3, sket.pnt.y[0]-3, 7, 7);

    for (var ipnt = 1; ipnt < npnt; ipnt++) {
        context.fillRect(sket.pnt.x[ipnt]-2, sket.pnt.y[ipnt]-2, 5, 5);
    }

    // display constraint labels associated with the Points
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (sket.pnt.lbl[ipnt] != "") {
            context.font         = "10px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            var width  = context.measureText(sket.pnt.lbl[ipnt]).width;
            var height = 8;

            context.fillStyle    = "yellow";
            context.fillRect(sket.pnt.x[ipnt]-width/2-2, sket.pnt.y[ipnt]-height/2-2, width+4, height+4);

            context.fillStyle    = "black";
            context.fillText(sket.pnt.lbl[ipnt], sket.pnt.x[ipnt], sket.pnt.y[ipnt]);
        }
    }

    // display constraint labels associated with the Segments
    for (var iseg = 0; iseg < nseg; iseg++) {
        if (sket.seg.lbl[iseg] != "") {
            context.font         = "12px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            var width  = context.measureText(sket.seg.lbl[iseg]).width;
            var height = 8;

            context.fillStyle    = "yellow";
            context.fillRect(sket.seg.xm[iseg]-width/2-2, sket.seg.ym[iseg]-height/2-2, width+4, height+4);

            context.fillStyle    = "black";
            context.fillText(sket.seg.lbl[iseg], sket.seg.xm[iseg], sket.seg.ym[iseg]);
        }
    }

    // display width and depth constraints
    for (var icon = 0; icon < ncon; icon++) {
        if (sket.con.type[icon] == "W" || sket.con.type[icon] == "D") {
            var ibeg = sket.con.index1[icon];
            var iend = sket.con.index2[icon];

            context.strokeStyle = "green";
            context.lineWidth   = 1;
            context.beginPath();
            context.moveTo(sket.pnt.x[ibeg], sket.pnt.y[ibeg]);
            if (sket.con.type[icon] == "W") {
                var xmid = (  sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 2;
                var ymid = (2*sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 3;

                context.lineTo(sket.pnt.x[ibeg], ymid);
                context.lineTo(sket.pnt.x[iend], ymid);
            } else {
                var xmid = (2*sket.pnt.x[ibeg] + sket.pnt.x[iend]) / 3;
                var ymid = (  sket.pnt.y[ibeg] + sket.pnt.y[iend]) / 2;

                context.lineTo(xmid, sket.pnt.y[ibeg]);
                context.lineTo(xmid, sket.pnt.y[iend]);
            }
            context.lineTo(sket.pnt.x[iend], sket.pnt.y[iend]);
            context.stroke();

            context.font         = "12px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            if (sket.con.type[icon] == "W") {
                var label = "W";
            } else {
                var label = "D";
            }
            var width  = context.measureText(label).width;
            var height = 8;

            context.fillStyle = "yellow";
            context.fillRect(xmid-width/2-2, ymid-height/2-2, width+4, height+4);

            context.fillStyle = "black";
            context.fillText(label, xmid, ymid);
        }
    }

    if (sket.suggest) {
        // suggested Sketch deletions
        if        (sket.suggest.substring(0,4) == "*del") {
            var suggestions = sket.suggest.split(";");

            for (ipnt = 0; ipnt < sket.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] == -1     ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "red";

                    context.clearRect(sket.pnt.x[ipnt]-width/2-2, sket.pnt.y[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sket.pnt.x[ipnt], sket.pnt.y[ipnt]-3*height/2);
                }
            }

            for (ipnt = 0; ipnt < sket.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] >  -1       ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "red";

                    context.clearRect(sket.seg.xm[ipnt]-width/2-2, sket.seg.ym[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sket.seg.xm[ipnt], sket.seg.ym[ipnt]-3*height/2);
                }
            }

            return;

        // suggested Sketch additions
        } else if (sket.suggest.substring(0,4) == "*add") {
            var suggestions = sket.suggest.split(";");

            for (ipnt = 0; ipnt < sket.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] == -1     ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "green";

                    context.clearRect(sket.pnt.x[ipnt]-width/2-2, sket.pnt.y[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sket.pnt.x[ipnt], sket.pnt.y[ipnt]-3*height/2);
                }
            }

            for (ipnt = 0; ipnt < sket.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] >  -1       ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "green";

                    context.clearRect(sket.seg.xm[ipnt]-width/2-2, sket.seg.ym[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sket.seg.xm[ipnt], sket.seg.ym[ipnt]-3*height/2);
                }
            }

            return;
        }
    }

    // if we are in mode 1, draw a line from the last Point to the cursor
    if (sket.mode == 1) {

        // nominal line to cursor is blue
        var color = "blue";

        // if horizontal, make it orange
        if (Math.abs(sket.pnt.y[npnt-1]-wv.cursorY) < sket.halo) {
            color = "orange";

        // otherwise see if it alignes with another point
        } else {
            for (var jpnt = 0; jpnt < npnt-1; jpnt++) {
                if (Math.abs(sket.pnt.y[jpnt]-wv.cursorY) < sket.halo) {
                    context.strokeStyle = "orange";
                    context.lineWidth   = 1;

                    context.beginPath();
                    context.moveTo(sket.pnt.x[jpnt], sket.pnt.y[jpnt]);
                    context.lineTo(wv.cursorX,  wv.cursorY );
                    context.stroke();
                    break;
                }
            }
        }

        // if vertical, make it orange
        if (Math.abs(sket.pnt.x[npnt-1]-wv.cursorX) < sket.halo) {
            color = "orange";

        // otherwise see if it alignes with another point
        } else {
            for (var jpnt = 0; jpnt < npnt-1; jpnt++) {
                if (Math.abs(sket.pnt.x[jpnt]-wv.cursorX) < sket.halo) {
                    context.strokeStyle = "orange";
                    context.lineWidth   = 1;

                    context.beginPath();
                    context.moveTo(sket.pnt.x[jpnt], sket.pnt.y[jpnt]);
                    context.lineTo(wv.cursorX,  wv.cursorY );
                    context.stroke();
                    break;
                }
            }
        }

        // if the cursor is close to the first Point, indicate that
        if (Math.abs(sket.pnt.x[0]-wv.cursorX) < sket.halo &&
            Math.abs(sket.pnt.y[0]-wv.cursorY) < sket.halo   ) {
            context.strokeStyle = "orange";
            context.lineWidth   = 5;
            context.beginPath();
            context.arc(sket.pnt.x[0], sket.pnt.y[0], 17, 0, 2*Math.PI);
            context.stroke();
        }

        // draw the line to the cursor
        context.strokeStyle = color;
        context.lineWidth   = 2;

        context.beginPath();
        context.moveTo(sket.pnt.x[npnt-1], sket.pnt.y[npnt-1]);
        context.lineTo(wv.cursorX,    wv.cursorY   );
        context.stroke();
    }

    // have the Build button show the current Sketcher status
    var button = document.getElementById("buildButton");
    button.style.backgroundColor =  "#FFFF3F";

    if        (sket.mode == 0) {
        button["innerHTML"] = "Initializing...";
    } else if (sket.mode == 1) {
        button["innerHTML"] = "Drawing...";
    } else if (sket.mode == 2) {
        button["innerHTML"] = "Setting R...";
//  } else if (sket.mode == 3) {
//      handled below
    } else if (sket.mode == 4) {
        button["innerHTML"] = "Setting W...";
    } else if (sket.mode == 5) {
        button["innerHTML"] = "Setting D...";
    } else if (sket.mode == 6) {
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;
        document.getElementById("sketchMenuBtn").style.backgroundColor = "#3FFF3F";
    } else if (sket.var.name.length != sket.con.type.length) {
        button["innerHTML"] = "Constraining...";
        button.style.backgroundColor = "#FFFF3F"
        document.getElementById("sketchMenuBtn").style.backgroundColor = null;
    } else {
        button["innerHTML"] = "Press to Solve";
        button.style.backgroundColor = "#3FFF3F";
    }

    // post informtion about current mode in blframe
    var skstat = document.getElementById("SketchStatus");

    var mesg = "ndof=" + sket.var.name.length + "   ncon=" + sket.con.type.length + "\n";
    if        (sket.mode == 1) {
        mesg += "Valid commands are:\n";
        mesg += "  'l'   add linseg\n";
        mesg += "  'c'   add cirarc\n";
        mesg += "  's'   add spline\n";
        mesg += "  'b'   add bezier\n";
        mesg += "  'z'   add zero-length segment";
        mesg += "  'o'   complete (open) sketch\n";
    } else if (sket.mode == 2) {
        mesg += "Hit any character to set curvature\n";
    } else if (sket.mode == 3) {
        mesg += "Valid constraints at points\n";
        mesg += "  'x' (fix x)     'y' (fix y)\n";
        mesg += "  'p' (perp)      't' (tangent)\n";
        mesg += "  'a' (angle)\n";
        mesg += "  'w' (width)     'd' (depth)\n";
        mesg += "Valid constraints on segments\n";
        mesg += "  'h' (horiz)     'v' (vertical)\n";
        mesg += "  'i' (incline)   'l' (length)\n";
        mesg += "Valid constraints on cirarcs\n";
        mesg += "  'r' (radius)    's' (sweep angle)\n";
        mesg += "Valid commands anywhere\n";
        mesg += "  '<' (delete)    '?' (info)\n";
    } else if (sket.mode == 4) {
        mesg += "Hit any character to set width\n";
    } else if (sket.mode == 5) {
        mesg += "Hit any character to set depth\n";
    }

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    skstat.replaceChild(pre, skstat.lastChild);
}


// draw a Bezier curve (with evaluations via Hoerner-like scheme)
//
function drawBezier(context, Xbezier, Ybezier)
{
    // alert("drawBezier()");

    var degree = Xbezier.length - 1;

    for (var j = 1; j <= 4*degree; j++) {
        var t  = j / (4.0 * degree);
        var t1 = 1.0 - t;

        var fact   = 1.0;
        var bicoef = 1;

        var xx = Xbezier[0] * t1;
        var yy = Ybezier[0] * t1;

        for (var i = 1; i < degree; i++) {
            fact    *= t;
            bicoef  *= (degree - i + 1) / i;

            xx = (xx + fact * bicoef * Xbezier[i]) * t1;
            yy = (yy + fact * bicoef * Ybezier[i]) * t1;
        }

        xx += fact * t * Xbezier[degree];
        yy += fact * t * Ybezier[degree];

        context.lineTo(xx, yy);
    }
}


// draw a Spline curve (with evaluations via Hoerner-like scheme)
//
function drawSpline(context, Xspline, Yspline)
{
    // alert("drawSpline()");

    var degree = Xspline.length - 1;

    // convert Spline points into cubic Bezier control points
    var Xbezier = [];
    var Ybezier = [];
    var gam     = [];

    var bet    = 4;
    Xbezier[1] = (6 * Xspline[1] - Xspline[0]) / bet;
    Ybezier[1] = (6 * Yspline[1] - Yspline[0]) / bet;

    for (var i = 1; i < degree-2; i++) {
        gam[i]       = 1 / bet;
        bet          =  4 - gam[i];
        Xbezier[i+1] = (6 * Xspline[i+1] - Xbezier[i]) / bet;
        Ybezier[i+1] = (6 * Yspline[i+1] - Ybezier[i]) / bet;
    }

    gam[degree-2]     = 1 / bet;
    bet               =  4 - gam[degree-2];
    Xbezier[degree-1] = (6 * Xspline[degree-1] - Xspline[degree] - Xbezier[degree-2]) / bet;
    Ybezier[degree-1] = (6 * Yspline[degree-1] - Yspline[degree] - Ybezier[degree-2]) / bet;

    for (i = degree-2; i > 0; i--) {
        Xbezier[i] -= gam[i] * Xbezier[i+1];
        Ybezier[i] -= gam[i] * Ybezier[i+1];
    }

    Xbezier[0] = Xspline[0];
    Ybezier[0] = Yspline[0];

    Xbezier[degree] = Xspline[degree];
    Ybezier[degree] = Yspline[degree];

    // draw the Splines
    for (i = 1; i <= degree; i++) {
        for (var j = 1; j < 5; j++) {
            var t  = 0.25 * j;
            var t1 = 1.0 - t;

            var xx = (t1 * t1 * t1) * (    Xspline[i-1]                 )
                   + (t1 * t1 * t ) * (2 * Xbezier[i-1] +     Xbezier[i])
                   + (t1 * t  * t ) * (    Xbezier[i-1] + 2 * Xbezier[i])
                   + (t  * t  * t ) * (                       Xspline[i]);
            var yy = (t1 * t1 * t1) * (    Yspline[i-1]                 )
                   + (t1 * t1 * t ) * (2 * Ybezier[i-1] +     Ybezier[i])
                   + (t1 * t  * t ) * (    Ybezier[i-1] + 2 * Ybezier[i])
                   + (t  * t  * t ) * (                       Yspline[i]);
            context.lineTo(xx,yy);
        }
    }
}


//
// get Point closest to the cursor
//
function getClosestSketchPoint()
{
    // alert("getClosestSketchPoint()");

    var npnt  = sket.pnt.x.length;
    var ibest = -1;
    var dbest = 999999;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var dtest = Math.abs(wv.cursorX - sket.pnt.x[ipnt])
                  + Math.abs(wv.cursorY - sket.pnt.y[ipnt]);
        if (dtest < dbest) {
            ibest = ipnt;
            dbest = dtest;
        }
    }

    return ibest;
}


//
// get Segment closest to the cursor
//
function getClosestSketchSegment()
{
    // alert("getClosestSketchSegment()");

    var nseg  = sket.seg.type.length;
    var ibest = -1;
    var dbest = 999999;

    for (var iseg = 0; iseg < nseg; iseg++) {
        var dtest = Math.abs(wv.cursorX - sket.seg.xm[iseg])
                  + Math.abs(wv.cursorY - sket.seg.ym[iseg]);
        if (dtest < dbest) {
            ibest = iseg;
            dbest = dtest;
        }
    }

    return ibest;
}


//
// get W or D constrint closest to the cursor
//
function getClosestSketchConstraint()
{
    // alert("getClosestSketchConstraint()");

    var ncon  = sket.con.type.length;
    var ibest = undefined;
    var dbest = 999999;

    for (var icon = 0; icon < ncon; icon++) {
        if        (sket.con.type[icon] == "W") {
            var xmid = (  sket.pnt.x[sket.con.index1[icon]]
                        + sket.pnt.x[sket.con.index2[icon]]) / 2;
            var ymid = (2*sket.pnt.y[sket.con.index1[icon]]
                        + sket.pnt.y[sket.con.index2[icon]]) / 3;
        } else if (sket.con.type[icon] == "D") {
            var xmid = (2*sket.pnt.x[sket.con.index1[icon]]
                        + sket.pnt.x[sket.con.index2[icon]]) / 3;
            var ymid = (  sket.pnt.y[sket.con.index1[icon]]
                        + sket.pnt.y[sket.con.index2[icon]]) / 2;
        } else {
            continue;
        }

        var dtest = Math.abs(wv.cursorX - xmid) + Math.abs(wv.cursorY - ymid);
        if (dtest < dbest) {
            ibest = icon;
            dbest = dtest;
        }
    }

    return ibest;
}


//
// constructor for a Tree
//
function Tree(doc, treeId) {
    // alert("in Tree(doc="+doc+", treeId="+treeId+")");

    // remember the document
    this.document = doc;
    this.treeId   = treeId;

    // arrays to hold the Nodes
    this.name    = new Array();
    this.tooltip = new Array();
    this.gprim   = new Array();
    this.click   = new Array();
    this.parent  = new Array();
    this.child   = new Array();
    this.next    = new Array();
    this.nprop   = new Array();
    this.opened  = new Array();

    this.prop1  = new Array();
    this.cbck1  = new Array();
    this.prop2  = new Array();
    this.cbck2  = new Array();
    this.prop3  = new Array();
    this.cbck3  = new Array();

    // initialize Node=0 (the root)
    this.name[  0]  = "**root**";
    this.tooltip[0] = "";
    this.gprim[ 0]  = "";
    this.click[ 0]  = null;
    this.parent[0]  = -1;
    this.child[ 0]  = -1;
    this.next[  0]  = -1;
    this.nprop[ 0]  =  0;
    this.prop1[ 0]  = "";
    this.cbck1[ 0]  = null;
    this.prop2[ 0]  = "";
    this.cbck2[ 0]  = null;
    this.prop3[ 0]  = "";
    this.cbck3[ 0]  = null;
    this.opened[0]  = +1;

    // add methods
    this.addNode  = TreeAddNode;
    this.expand   = TreeExpand;
    this.contract = TreeContract;
    this.prop     = TreeProp;
    this.clear    = TreeClear;
    this.build    = TreeBuild;
    this.update   = TreeUpdate;
}


//
// add a Node to the Tree
//
function TreeAddNode(iparent, name, tooltip, gprim, click,
                     prop1, cbck1,
                     prop2, cbck2,
                     prop3, cbck3)
{
    // alert("in TreeAddNode(iparent="+iparent+", name="+name+", tooltip="+tooltip+", gprim="+gprim+
    //                       ", click="+click+
    //                       ", prop1="+prop1+", cbck1="+cbck1+
    //                       ", prop2="+prop2+", cbck2="+cbck2+
    //                       ", prop3="+prop3+", cbck3="+cbck3+")");

    // validate the input
    if (iparent < 0 || iparent >= this.name.length) {
        alert("iparent="+iparent+" is out of range");
        return;
    }

    // find the next Node index
    var inode = this.name.length;

    // store this Node's values
    this.name[   inode] = name;
    this.tooltip[inode] = tooltip;
    this.gprim[  inode] = gprim;
    this.click[  inode] = click;
    this.parent[ inode] = iparent;
    this.child[  inode] = -1;
    this.next[   inode] = -1;
    this.nprop[  inode] =  0;
    this.opened[ inode] =  0;

    // store the properties
    if (prop1 !== undefined) {
        this.nprop[inode] = 1;
        this.prop1[inode] = prop1;
        this.cbck1[inode] = cbck1;
    }

    if (prop2 !== undefined) {
        this.nprop[inode] = 2;
        this.prop2[inode] = prop2;
        this.cbck2[inode] = cbck2;
    }

    if (prop3 !== undefined) {
        this.nprop[inode] = 3;
        this.prop3[inode] = prop3;
        this.cbck3[inode] = cbck3;
    }

    // if the parent does not have a child, link this
    //    new Node to the parent
    if (this.child[iparent] < 0) {
        this.child[iparent] = inode;

    // otherwise link this Node to the last parent's child
    } else {
        var jnode = this.child[iparent];
        while (this.next[jnode] >= 0) {
            jnode = this.next[jnode];
        }

        this.next[jnode] = inode;
    }
}


//
// build the Tree (ie, create the html table from the Nodes)
//
function TreeBuild() {
    // alert("in TreeBuild()");

    var doc = this.document;

    // if the table already exists, delete it and all its children (3 levels)
    var thisTable = doc.getElementById(this.treeId);
    if (thisTable) {
        var child1 = thisTable.lastChild;
        while (child1) {
            var child2 = child1.lastChild;
            while (child2) {
                var child3 = child2.lastChild;
                while (child3) {
                    child2.removeChild(child3);
                    child3 = child2.lastChild;
                }
                child1.removeChild(child2);
                child2 = child1.lastChild;
            }
            thisTable.removeChild(child1);
            child1 = thisTable.lastChild;
        }
        thisTable.parentNode.removeChild(thisTable);
    }

    // build the new table
    var newTable = doc.createElement("table");
    newTable.setAttribute("id", this.treeId);
    doc.getElementById("treefrm").appendChild(newTable);

    // traverse the Nodes using depth-first search
    var inode = 1;
    while (inode > 0) {

        // table row "node"+inode
        var newTR = doc.createElement("TR");
        newTR.setAttribute("id", "node"+inode);
        newTable.appendChild(newTR);

        // table data "node"+inode+"col1"
        var newTDcol1 = doc.createElement("TD");
        newTDcol1.setAttribute("id", "node"+inode+"col1");
        newTR.appendChild(newTDcol1);

        var newTexta = doc.createTextNode("");
        newTDcol1.appendChild(newTexta);

        // table data "node"+inode+"col2"
        var newTDcol2 = doc.createElement("TD");
        newTDcol2.setAttribute("id", "node"+inode+"col2");
        if (this.click[inode] != null) {
            newTDcol2.className = "fakelinkcmenu";
            if (this.tooltip[inode].length > 0) {
                newTDcol2.title = this.tooltip[inode];
            }
        }
        newTR.appendChild(newTDcol2);

        var newTextb = doc.createTextNode(this.name[inode]);
        newTDcol2.appendChild(newTextb);

        var name = this.name[inode].replace(/\u00a0/g, "");

        var ibrch = 0;
        for (var jbrch = 0; jbrch < brch.length; jbrch++) {
            if (brch[jbrch].name == name) {
                if (brch[jbrch].ileft == -2) {
                    newTDcol2.className = "errorTD";
                }
                break;
            }
        }

        // table data "node"+inode+"col3"
        if (this.nprop[inode] > 0) {
            var newTDcol3 = doc.createElement("TD");
            newTDcol3.setAttribute("id", "node"+inode+"col3");
            if (this.cbck1[inode] != "") {
                newTDcol3.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol3);

            if (this.nprop[inode] == 1) {
                newTDcol3.setAttribute("colspan", "3");
            }

            var newTextc = doc.createTextNode(this.prop1[inode]);
            newTDcol3.appendChild(newTextc);
        }

        // table data "node:+inode+"col4"
        if (this.nprop[inode] > 1) {
            var newTDcol4 = doc.createElement("TD");
            newTDcol4.setAttribute("id", "node"+inode+"col4");
            if (this.cbck2[inode] != "") {
                newTDcol4.className = "fakelinkoff";
            }
            newTR.appendChild(newTDcol4);

            if (this.nprop[inode] == 2) {
                newTDcol4.setAttribute("colspan", "2");
            }

            var newTextd = doc.createTextNode(this.prop2[inode]);
            newTDcol4.appendChild(newTextd);
        }

        // table data "node:+inode+"col5"
        if (this.nprop[inode] > 2) {
            var newTDcol5 = doc.createElement("TD");
            newTDcol5.setAttribute("id", "node"+inode+"col5");
            if (this.cbck3[inode] != "") {
                newTDcol5.className = "fakelinkoff";
            }
            newTR.appendChild(newTDcol5);

            var newTextd = doc.createTextNode(this.prop3[inode]);
            newTDcol5.appendChild(newTextd);
        }

        // go to next row
        if        (this.child[inode] >= 0) {
            inode = this.child[inode];
        } else if (this.next[inode] >= 0) {
            inode = this.next[inode];
        } else {
            while (inode > 0) {
                inode = this.parent[inode];
                if (this.parent[inode] == 0) {
                    newTR = doc.createElement("TR");
                    newTR.setAttribute("height", "10px");
                    newTable.appendChild(newTR);
                }
                if (this.next[inode] >= 0) {
                    inode = this.next[inode];
                    break;
                }
            }
        }
    }

    this.update();
}


//
// clear the Tree
//
function TreeClear() {
    // alert("in TreeClear()");

    // remove all but the first Node
    this.name.splice(   1);
    this.tooltip.splice(1);
    this.gprim.splice(  1);
    this.click.splice(  1);
    this.parent.splice( 1);
    this.child.splice(  1);
    this.next.splice(   1);
    this.nprop.splice(  1);
    this.opened.splice( 1);

    this.prop1.splice(1);
    this.cbck1.splice(1);
    this.prop2.splice(1);
    this.cbck2.splice(1);
    this.prop3.splice(1);
    this.cbck3.splice(1);

    // reset the root Node
    this.parent[0] = -1;
    this.child[ 0] = -1;
    this.next[  0] = -1;
}


//
// expand a Node in the Tree
//
function TreeContract(inode) {
    // alert("in TreeContract(inode="+inode+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // contract inode
    this.opened[inode] = 0;

    // contract all descendents of inode
    for (var jnode = 1; jnode < this.parent.length; jnode++) {
        var iparent = this.parent[jnode];
        while (iparent > 0) {
            if (iparent == inode) {
                this.opened[jnode] = 0;
                break;
            }

            iparent = this.parent[iparent];
        }
    }

    // update the display
    this.update();
}


//
// expand a Node in the Tree
//
function TreeExpand(inode) {
    // alert("in TreeExpand(inode="+inode+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // expand inode
    this.opened[inode] = 1;

    // update the display
    this.update();
}


//
// change a property of a Node
//
function TreeProp(inode, iprop, onoff) {
    // alert("in TreeProp(inode="+inode+", iprop="+iprop+", onoff="+onoff+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    } else if (onoff != "on" && onoff != "off") {
        alert("onoff="+onoff+" is not 'on' or 'off'");
        return;
    }

    if (this != myTree) {
        alert("this="+this+"   myTree="+myTree);
    }

    var thisNode = "";

    // set the property for inode
    if        (iprop == 1 && this.prop1[inode] == "Viz") {
        thisNode = this.document.getElementById("node"+inode+"col3");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ON;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ON;
            }
        }
    } else if (iprop == 1) {
    } else if (iprop == 2 && this.prop2[inode] == "Grd") {
        thisNode = this.document.getElementById("node"+inode+"col4");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.LINES;
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.POINTS;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.LINES;
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.POINTS;
            }
        }
    } else if (iprop == 2) {
    } else if (iprop == 3 && this.prop3[inode] == "Trn") {
        thisNode = this.document.getElementById("node"+inode+"col5");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.TRANSPARENT;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.TRANSPARENT;
            }
        }
    } else if (iprop == 3 && this.prop3[inode] == "Ori") {
        thisNode = this.document.getElementById("node"+inode+"col5");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ORIENTATION;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ORIENTATION;
            }
        }
    } else if (iprop ==3) {
    } else {
        alert("iprop="+iprop+" is not 1, 2, or 3");
        return;
    }

    // update fakelinks in TreeWindow (needed when .attrs do not exist)
    if (thisNode != "") {
        if (onoff == "on") {
            thisNode.setAttribute("class", "fakelinkon");
            thisNode.title = "Toggle Ori off";
        } else {
            thisNode.setAttribute("class", "fakelinkoff");
            thisNode.title = "Toggle Ori on";
        }
    }

    // set property for inode's children
    for (var jnode = inode+1; jnode < this.parent.length; jnode++) {
        if (this.parent[jnode] == inode) {
            this.prop(jnode, iprop, onoff);
        }
    }

    wv.sceneUpd = 1;
}


//
// update the Tree (after build/expension/contraction/property-set)
//
function TreeUpdate() {
    // alert("in TreeUpdate()");

    var doc = this.document;

    // traverse the Nodes using depth-first search
    for (var inode = 1; inode < this.opened.length; inode++) {
        var element = doc.getElementById("node"+inode);

        // unhide the row
        element.style.display = "table-row";

        // hide the row if one of its parents has .opened=0
        var jnode = this.parent[inode];
        while (jnode != 0) {
            if (this.opened[jnode] == 0) {
                element.style.display = "none";
                break;
            }

            jnode = this.parent[jnode];
        }

        // if the current Node has children, set up appropriate event handler to expand/collapse
        if (this.child[inode] > 0) {
            if (this.opened[inode] == 0) {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "+";
                myElem.title   = "Expand";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.expand(thisNode);
                };

            } else {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "-";
                myElem.title   = "Collapse";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.contract(thisNode);
                };
            }
        }

        if (this.click[inode] !== null) {
            var myElem = doc.getElementById("node"+inode+"col2");
            myElem.onclick = this.click[inode];
        }

        // set the class of the properties
        if (this.nprop[inode] >= 1) {
            var myElem = doc.getElementById("node"+inode+"col3");
            myElem.onclick = this.cbck1[inode];

            if (this.prop1[inode] == "Viz") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Viz on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Viz off";
                    }
                }
            }
        }

        if (this.nprop[inode] >= 2) {
            var myElem = doc.getElementById("node"+inode+"col4");
            myElem.onclick = this.cbck2[inode];

            if (this.prop2[inode] == "Grd") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Grd on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Grd off";
                    }
                }
            }
        }

        if (this.nprop[inode] >= 3) {
            var myElem = doc.getElementById("node"+inode+"col5");
            myElem.onclick = this.cbck3[inode];

            if (this.prop3[inode] == "Trn") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Trn on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Trn off";
                    }
                }
            } else if (this.prop3[inode] == "Ori") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Ori on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Ori off";
                    }
                }
            }
        }
    }
}


//
// callback when "onresize" event occurs (called by ESP.html)
//
// resize the frames (with special handling to width of tlframe and height of brframe)
//
function resizeFrames() {
    // alert("resizeFrames()");

    var scrollSize = 24;

    // get the size of the client (minus amount to account for scrollbars)
    var body = document.getElementById("mainBody");
    var bodyWidth  = body.clientWidth  - scrollSize;
    var bodyHeight = body.clientHeight - scrollSize;

    // get the elements associated with the frames and the canvas
    var topleft = document.getElementById("tlframe");
    var butnfrm = document.getElementById("butnfrm");
    var treefrm = document.getElementById("treefrm");
    var toprite = document.getElementById("trframe");
    var botleft = document.getElementById("blframe");
    var botrite = document.getElementById("brframe");
    var canvas  = document.getElementById(wv.canvasID);
    var sketch  = document.getElementById("sketcher");

    // compute and set the widths of the frames
    //    (do not make tlframe larger than 250px)
    var leftWidth = Math.round(0.25 * bodyWidth);
    if (leftWidth > 250)   leftWidth = 250;
    var riteWidth = bodyWidth - leftWidth;
    var canvWidth = riteWidth - scrollSize;

    topleft.style.width = leftWidth+"px";
    butnfrm.style.width = leftWidth+"px";
    treefrm.style.width = leftWidth+"px";
    toprite.style.width = riteWidth+"px";
    botleft.style.width = leftWidth+"px";
    botrite.style.width = riteWidth+"px";
    canvas.style.width  = canvWidth+"px";
    canvas.width        = canvWidth;
    sketch.style.width  = canvWidth+"px";
    sketch.width        = canvWidth;

    // compute and set the heights of the frames
    //    (do not make brframe larger than 200px)
    var botmHeight = Math.round(0.20 * bodyHeight);
    if (botmHeight > 200)   botmHeight = 200;
    var  topHeight = bodyHeight - botmHeight;
    var canvHeight =  topHeight - scrollSize - 5;
    var keyHeight  = botmHeight - 25;

    topleft.style.height =  topHeight+"px";
    treefrm.style.height = (topHeight-90)+"px";
    toprite.style.height =  topHeight+"px";
    botleft.style.height = botmHeight+"px";
    botrite.style.height = botmHeight+"px";
    canvas.style.height  = canvHeight+"px";
    canvas.height        = canvHeight;
    sketch.style.height  = canvHeight+"px";
    sketch.height        = canvHeight;

    // set up canvas associated with key
    if (wv.canvasKY !== undefined) {
        var canf = document.getElementById(wv.canvasKY);

        var keyfWidth     = leftWidth - 20;
        canf.style.width  = keyfWidth + "px";
        canf.width        = keyfWidth;

        var keyfHeight    = botmHeight - 25;
        canf.style.height = keyfHeight + "px";
        canf.height       = keyfHeight;

        // force a key redraw
        wv.drawKey = 1;
    }
}


//
// change mode for trframe
//
function changeMode(newMode) {
    // alert("in changeMode(newMode="+newMode+")");

    var webViewer      = document.getElementById("WebViewer");
    var addBrchForm    = document.getElementById("addBrchForm");
    var editBrchForm   = document.getElementById("editBrchForm");
    var addBrchHeader  = document.getElementById("addBrchHeader");
    var editBrchHeader = document.getElementById("editBrchHeader");
    var editPmtrForm   = document.getElementById("editPmtrForm");
    var addPmtrHeader  = document.getElementById("addPmtrHeader");
    var editPmtrHeader = document.getElementById("editPmtrHeader");
    var editCsmForm    = document.getElementById("editCsmForm");
    var sketcherForm   = document.getElementById("sketcherForm");

    var wvKey          = document.getElementById("WVkey");
    var sketchStatus   = document.getElementById("SketchStatus");
    var ESPlogo        = document.getElementById("ESPlogo");

    if (newMode == wv.curMode) {
        return;
    } else if (newMode < 0) {
        // used to cause buttons such as "cmdSave" not to be active
        wv.curMode = newMode;

        var button = document.getElementById("stepThruBtn");
        button["innerHTML"] = "StepThru";

        return;
    } else if (newMode == 0) {
        webViewer.hidden    = false;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        wv.curMode   = 0;
        wv.curPmtr   = -1;
        wv.curBrch   = -1;
        wv.afterBrch = -1;
        wv.menuEvent = undefined;
    } else if (newMode == 1) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = false;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        // unselect all items
        var elements = document.getElementById("addBrchForm").elements;
        for (var ielem = 0; ielem < elements.length; ielem++) {
            elements[ielem].checked = false;
        }

        wv.curMode = 1;
    } else if (newMode == 2) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = false;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        addBrchHeader.hidden  = false;
        editBrchHeader.hidden = true;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curMode = 2;
    } else if (newMode == 3) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = false;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        addBrchHeader.hidden  = true;
        editBrchHeader.hidden = false;

        showBrchArgs();

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curMode = 3;
    } else if (newMode == 4) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = false;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        addPmtrHeader.hidden  = false;
        editPmtrHeader.hidden = true;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curMode = 4;
    } else if (newMode == 5) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = false;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        addPmtrHeader.hidden  = true;
        editPmtrHeader.hidden = false;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curMode = 5;
    } else if (newMode == 6) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = false;
        sketcherForm.hidden = true;
        wvKey.hidden        = true;
        sketchStatus.hidden = true;
        ESPlogo.hidden      = false;

        wv.curMode = 6;
    } else if (newMode == 7) {
        webViewer.hidden    = true;
        addBrchForm.hidden  = true;
        editBrchForm.hidden = true;
        editPmtrForm.hidden = true;
        editCsmForm.hidden  = true;
        sketcherForm.hidden = false;
        wvKey.hidden        = true;
        sketchStatus.hidden = false;
        ESPlogo.hidden      = true;

        wv.curMode = 7;
    } else {
        alert("Bad new mode = "+newMode);
    }

    if        (wv.curMode == 0) {
        document.getElementById("sketchMenuBtn").hidden = true;
    } else if (wv.curMode == 7) {
        document.getElementById("sketchMenuBtn").hidden = false;
    }
}


//
// rebuild the Tree Window
//
function rebuildTreeWindow() {
    // alert("in rebuildTreeWindow()");

    // if there was a previous Tree, keep track of whether or not
    //    the Parameters, Branches, and Display was open
    var pmtr1Open = 0;
    var pmtr2Open = 0;
    var brchsOpen = 0;

    if (myTree.opened.length > 4) {
        pmtr1Open = myTree.opened[1];
        pmtr2Open = myTree.opened[2];
        brchsOpen = myTree.opened[3];
    }

    // clear previous Nodes from the Tree
    myTree.clear();

    // put the group headers into the Tree
    myTree.addNode(0, "Design Parameters", "Add a Parameter",   "", addPmtr);
    myTree.addNode(0, "Local Variables",   "",                  "", null   );
    myTree.addNode(0, "Branches",          "Add Branch at end", "", addBrch);
    if (wv.curStep == 0) {
        myTree.addNode(0, "Display",           "",                  "", null,
                       "Viz", toggleViz,
                       "Grd", toggleGrd);
    } else {
        myTree.addNode(0, "Display");
    }

    // put the Design Parameters into the Tree
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        var name  = "\u00a0\u00a0"+pmtr[ipmtr].name;
        var type  =                pmtr[ipmtr].type;
        var nrow  =                pmtr[ipmtr].nrow;
        var ncol  =                pmtr[ipmtr].ncol;
        var value =                pmtr[ipmtr].value[0];

        if (nrow > 1 || ncol > 1) {
            value = "["+nrow+"x"+ncol+"]";
        }

        if (type == 500 && pmtr[ipmtr].dot[0] != 0) {
            myTree.addNode(1, "^"+name, "Edit Parameter", "", editPmtr,
                           ""+value, "");
        } else if (type == 500) {
            myTree.addNode(1,     name, "Edit Parameter", "", editPmtr,
                           ""+value, "");
        }
    }

    wv.pmtrStat = -2;

    // open the Design Parameters (if they were open before the Tree was rebuilt)
    if (pmtr1Open == 1) {
        myTree.opened[1] = 1;
    }

    // put the Local Variables into the Tree
    for (ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        var name  = "\u00a0\u00a0"+pmtr[ipmtr].name;
        var type  =                pmtr[ipmtr].type;
        var nrow  =                pmtr[ipmtr].nrow;
        var ncol  =                pmtr[ipmtr].ncol;
        var value =                pmtr[ipmtr].value[0];

        if (nrow > 1 || ncol > 1) {
            value = "["+nrow+"x"+ncol+"]";
        }

        if (type !== 500) {
            myTree.addNode(2, name, "", "", null,
                           ""+value, "");
            var inode = myTree.name.length - 1;

            var indx = 0;
            if (nrow > 1 || ncol > 1) {
                for (var irow = 0; irow < nrow; irow++) {
                    for (var icol = 0; icol < ncol; icol++) {
                        name  = "\u00a0\u00a0\u00a0\u00a0["+(irow+1)+","+(icol+1)+"]";
                        value = pmtr[ipmtr].value[indx++];

                        myTree.addNode(inode, name, "", "", null,
                                       "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0"+value, "");
                    }
                }
            }
        }
    }

    wv.pmtrStat = -2;

    // open the Local Variables (if they were open before the Tree was rebuilt)
    if (pmtr2Open == 1) {
        myTree.opened[2] = 1;
    }

    // put the Branches into the Tree
    var parents = [3];
    for (var ibrch = 0; ibrch < brch.length; ibrch++) {
        var name  = "\u00a0\u00a0";
        for (var indent = 0; indent < brch[ibrch].indent; indent++) {
            name = name+">";
        }
        name = name+brch[ibrch].name;

        if (ibrch == 0) {
        } else if (brch[ibrch].indent > brch[ibrch-1].indent) {
            parents.push(myTree.name.length-1);
        } else if (brch[ibrch].indent < brch[ibrch-1].indent) {
            parents.pop();
        }

        var type  =      brch[ibrch].type;
        var actv;
        if (ibrch >= wv.builtTo) {
            actv = "skipped";
        } else if (brch[ibrch].actv == 301) {
            actv = "suppressed";
        } else if (brch[ibrch].actv == 302) {
            actv = "inactive";
        } else if (brch[ibrch].actv == 303) {
            actv = "deferred";
        } else {
            actv = "";
        }

        myTree.addNode(parents[parents.length-1], name, "Edit/del/add-after Branch", "", editBrch,
                       type, "",
                       actv, "");
    }
    parents = undefined;

    wv.brchStat = -2;

    // open the Branches (if they were open before the Tree was rebuilt)
    if (brchsOpen == 1) {
        myTree.opened[3] = 1;
    }

    // put the Display attributes into the Tree
    var patchesNode = -1;    // tree Node that will contain the Patches
    for (var gprim in wv.sceneGraph) {
        if (wv.curStep == 1) {
            myTree.addNode(4, "\u00a0\u00a0CancelStepThru", "Cancel StepThru mode", null, cancelStepThru);
            break;
        }

        // parse the name
        var matches = gprim.split(" ");

        var iface = -1;
        var iedge = -1;
        var csys  = -1;

        if        (matches[0] == "Axes") {
            myTree.addNode(4, "\u00a0\u00a0Axes", "", gprim, null,
                           "Viz", toggleViz);
            myTree.addNode(4, "\u00a0\u00a0DisplayType", "Modify display type",     null, modifyDisplayType);
            myTree.addNode(4, "\u00a0\u00a0DisplayFilter", "Modify display filter", null, modifyDisplayFilter);
            continue;

        // processing for a Patch: "Patch m @I=n"
        } else if (matches.length == 3 && matches[0] == "Patch" &&
                   (matches[2].includes("@I=") ||
                    matches[2].includes("@J=") ||
                    matches[2].includes("@K=")   )                 ) {
            if (patchesNode == -1) {
                myTree.addNode(4, "\u00a0\u00a0Patches", "", gprim, null,
                              "Viz", toggleViz);
                patchesNode = myTree.name.length - 1;
            }
            myTree.addNode(patchesNode, "\u00a0\u00a0\u00a0\u00a0"+gprim, "", gprim, null,
                           "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing when Body is explicitly named: "Bodyname"
        } else if (matches.length == 1) {
            var bodyName = matches[0];
            var iface    = -1;
            var iedge    = -1;
            var csys     = -1;

        // processing when Body is not explicitly named: "Body m"
        } else if (matches.length == 2 && matches[0] == "Body") {
            var bodyName = matches[0] + " " + matches[1];
            var iface    = -1;
            var iedge    = -1;
            var csys     = -1;

        // processing when Body is explicitly named: "Bodyname Edge m"
        } else if (matches.length == 3) {
            var bodyName = matches[0];
            if        (matches[1] == "Face") {
                var iface = matches[2];
            } else if (matches[1] == "Edge") {
                var iedge = matches[2];
            } else if (matches[1] == "Csys") {
                var csys = matches[2];
            }

        // processing when Body is not explicitly named: "Body m Edge n"
        } else if (matches.length == 4 && matches[0] == "Body") {
            var bodyName = matches[0] + " " + matches[1];
            if        (matches[2] == "Face") {
                var iface = matches[3];
            } else if (matches[2] == "Edge") {
                var iedge = matches[3];
            } else if (matches[2] == "Csys") {
                var csys  = matches[3];
            }
        }

        // determine if Body does not exists
        var kbody = -1;
        for (var jnode = 1; jnode < myTree.name.length; jnode++) {
            if (myTree.name[jnode] == "\u00a0\u00a0"+bodyName) {
                kbody = jnode;
            }
        }

        // if Body does not exist, create it and its Face, Edge, and Csystem lists
        //    subnodes now
        var kface, kedge, kcsys;
        if (kbody < 0) {
            myTree.addNode(4, "\u00a0\u00a0"+bodyName, "Show Body Attributes", "", showBodyAttrs,
                           "Viz", toggleViz,
                           "Grd", toggleGrd);
            kbody = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Faces", "", "", null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Trn", toggleTrn);
            kface = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Edges", "", "", null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Ori", toggleOri);
            kedge = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Csystems", "", "", null,
                           "Viz", toggleViz);
            kcsys = myTree.name.length - 1;

        // otherwise, get pointers to the face-group and edge-group Nodes
        } else {
            kface = myTree.child[kbody];
            kedge = kface + 1;
            kcsys = kedge + 1;
        }

        // make the Tree Node
        if        (iface >= 0) {
            myTree.addNode(kface, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0face "+iface, "", gprim, null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Trn", toggleTrn);
        } else if (iedge >= 0) {
            myTree.addNode(kedge, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0edge "+iedge, "", gprim, null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Ori", toggleOri);
        } else if (csys  >= 0) {
            myTree.addNode(kcsys, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0"+csys, "", gprim, null,
                           "Viz", toggleViz);
        }
    }

    // open the Display (by default)
    myTree.opened[4] = 1;

    // mark that we have (re-)built the Tree
    wv.sgUpdate = 0;

    // convert the abstract Tree Nodes into an HTML table
    myTree.build();
}


//
// post a message into the brframe
//
function postMessage(mesg) {
    // alert("in postMessage(mesg="+mesg+")");

    if (wv.debugUI) {
        console.log("postMessage: "+mesg.substring(0,40));
    }

    var botm = document.getElementById("brframe");

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);
    botm.insertBefore(pre, botm.lastChild);

    pre.scrollIntoView();
}


//
// load info into editBrchForm
//
function setupEditBrchForm() {
    // alert("setupEditBrchForm()");

    var ibrch = wv.curBrch;
    var name  = brch[ibrch].name;
    var type  = brch[ibrch].type;
    var level = brch[ibrch].level;

    // return if within a UDC
    if (level > 0) {
        return 1;
    }

    var editBrchForm = document.getElementById("editBrchForm");

    // turn all arguments off (by default)
    document.getElementById("argName1").parentNode.style.display = "none";
    document.getElementById("argName2").parentNode.style.display = "none";
    document.getElementById("argName3").parentNode.style.display = "none";
    document.getElementById("argName4").parentNode.style.display = "none";
    document.getElementById("argName5").parentNode.style.display = "none";
    document.getElementById("argName6").parentNode.style.display = "none";
    document.getElementById("argName7").parentNode.style.display = "none";
    document.getElementById("argName8").parentNode.style.display = "none";
    document.getElementById("argName9").parentNode.style.display = "none";

    // turn all attributes/csystems off (by default)
    document.getElementById("attrName1").parentNode.style.display = "none";
    document.getElementById("attrName2").parentNode.style.display = "none";
    document.getElementById("attrName3").parentNode.style.display = "none";
    document.getElementById("attrName4").parentNode.style.display = "none";
    document.getElementById("attrName5").parentNode.style.display = "none";
    document.getElementById("attrName6").parentNode.style.display = "none";
    document.getElementById("attrName7").parentNode.style.display = "none";
    document.getElementById("attrName8").parentNode.style.display = "none";
    document.getElementById("attrName9").parentNode.style.display = "none";

    // start by looking at arguments
    document.getElementById("editBrchArgs" ).hidden = false;
    document.getElementById("editBrchAttrs").hidden = true;

    // fill in the Branch type and name
    document.getElementById("brchType").firstChild["data"] = type;
    editBrchForm.            brchName.value                = name;

    // by default, the "Enter Sketcher" button is turned off
    document.getElementById("EnterSketcher").style.display = "none";

    // set up arguments based upon the type
    var argList;
    var defValue;
    var suppress = 0;   // set to 1 if type can be suppressed

    if        (type == "applycsys") {
        argList  = ["$csysName", "ibody"];
        defValue = ["",          "0"    ];
    } else if (type == "arc") {
        argList  = ["xend", "yend", "zend", "dist", "$plane"];
        defValue = ["",     "",     "",     "",     "xy"    ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "assert") {
        argList  = ["arg1", "arg2", "toler", "verify"];
        defValue = ["",     "",     "0",     "0"     ];
    } else if (type == "bezier") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "blend") {
        argList  = ["begList", "endList", "reorder", "oneFace"];
        defValue = ["0",       "0",       "0",       "0"      ];
        suppress = 1;
    } else if (type == "box") {
        argList  = ["xmin", "ymin", "zmin", "dx", "dy", "dz"];
        defValue = ["",     "",     "",     "",   "",   ""  ];
        suppress = 1;
    } else if (type == "catbeg") {
        argList  = ["errCode"];
        defValue = [""       ];
    } else if (type == "catend") {
        argList  = [];
        defValue = [];
    } else if (type == "chamfer") {
        argList  = ["radius", "edgeList"];
        defValue = ["",       "0"       ];
        suppress = 1;
    } else if (type == "cirarc") {
        argList  = ["xon", "yon", "zon", "xend", "yend", "zend"];
        defValue = ["",    "",    "",    "",     "",     ""    ];
    } else if (type == "combine") {
        argList  = [];
        defValue = [];
    } else if (type == "cone") {
        argList  = ["xvrtx", "yvrtx", "zvrtx", "xbase", "ybase", "zbase", "radius"];
        defValue = ["",      "",      "",      "",      "",      "",      ""      ];
        suppress = 1;
    } else if (type == "connect") {
        argList  = ["faceList1", "faceList2"];
        defValue = ["",          ""         ];
    } else if (type == "cylinder") {
        argList  = ["xbeg", "ybeg", "zbeg", "xend", "yend", "zend", "radius"];
        defValue = ["",     "",     "",     "",     "",     "",     ""      ];
        suppress = 1;
    } else if (type == "dimension") {
        argList  = ["$pmtrName", "nrow", "ncol"];
        defValue = ["",          "",     ""    ];
    } else if (type == "dump") {
        argList  = ["$filename", "remove", "tomark"];
        defValue = ["",          "0",      "0"     ];
    } else if (type == "else") {
        argList  = [];
        defValue = [];
    } else if (type == "elseif") {
        argList  = ["val1", "$op1", "val2", "$op2", "val3", "$op3", "val4"];
        defValue = ["",     "",     "",     "and",  "0",    "eq",   "0"   ];
//  } else if (type == "end") {
    } else if (type == "endif") {
        argList  = [];
        defValue = [];
    } else if (type == "extract") {
        argList  = ["index"];
        defValue = [""     ];
    } else if (type == "extrude") {
        argList  = ["dx", "dy", "dz"];
        defValue = ["",   "",   ""  ];
        suppress = 1;
    } else if (type == "fillet") {
        argList  = ["radius", "edgeList"];
        defValue = ["",       "0"       ];
        suppress = 1;
    } else if (type == "getattr") {
        argList  = ["$pmtrName", "attrID"];
        defValue = ["",          ""      ];
    } else if (type == "group") {
        argList  = [];
        defValue = [];
    } else if (type == "hollow") {
        argList  = ["thick",  "entList", "listStyle"];
        defValue = ["0",      "0"      , "0"        ];
        suppress = 1;
    } else if (type == "ifthen") {
        argList  = ["val1", "$op1", "val2", "$op2", "val3", "$op3", "val4"];
        defValue = ["",     "",     "",     "and",  "0",    "eq",   "0"   ];
    } else if (type == "import") {
        argList  = ["$filename", "bodynumber"];
        defValue = ["",          "1"         ];
        suppress = 1;
//  } else if (type == "interface") {
    } else if (type == "intersect") {
        argList  = ["$order", "index", "maxtol"];
        defValue = ["none",   "1",     "0"     ];
    } else if (type == "join") {
        argList  = ["toler"];
        defValue = ["0"    ];
    } else if (type == "linseg") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "loft") {
        alert("Consider using 'rule' or 'blend'");
        argList  = ["smooth"];
        defValue = [""      ];
        suppress = 1;
    } else if (type == "macbeg") {
        alert("'macbeg' is deprecated");
        argList  = ["imacro"];
        defValue = [""      ];
    } else if (type == "macend") {
        alert("'macend' is deprecated");
        argList  = [];
        defValue = [];
    } else if (type == "mark") {
        argList  = [];
        defValue = [];
    } else if (type == "mirror") {
        argList  = ["nx", "ny", "nz", "dist"];
        defValue = ["",   "",   "",   "0"   ];
        suppress = 1;
    } else if (type == "patbeg") {
        argList  = ["$pmtrName", "ncopy"];
        defValue = ["",          ""     ];
    } else if (type == "patbreak") {
        argList  = ["expr"];
        defValue = [""    ];
    } else if (type == "patend") {
        argList  = [];
        defValue = [];
    } else if (type == "point") {
        argList  = ["xloc", "yloc", "zloc"];
        defValue = ["",     "",     ""    ];
        suppress = 1;
    } else if (type == "project") {
        argList  = ["x", "y", "z", "dx", "dy", "dz", "useEdges"];
        defValue = ["",  "",  "",  "",   "",   "",   "0"       ];
    } else if (type == "recall") {
        aleft("'recall' is deprecated");
        argList  = ["imicro"];
        defValue = [""      ];
    } else if (type == "reorder") {
        argList  = ["ishift", "iflip"];
        defValue = ["",       "0"    ];
        suppress = 1;
    } else if (type == "restore") {
        argList  = ["$name", "index"];
        defValue = ["",      "0"    ];
    } else if (type == "revolve") {
        argList  = ["xorig", "yorig", "zorig", "dxaxis", "dyaxis", "dzaxis", "angDeg"];
        defValue = ["",      "",      "",      "",       "",       "",       ""      ];
        suppress = 1;
    } else if (type == "rotatex") {
        argList  = ["angDeg", "yaxis", "zaxis"];
        defValue = ["",       "",      ""     ];
        suppress = 1;
    } else if (type == "rotatey") {
        argList  = ["angDeg", "zaxis", "xaxis"];
        defValue = ["",       "",      ""     ];
        suppress = 1;
    } else if (type == "rotatez") {
        argList  = ["angDeg", "xaxis", "yaxis"];
        defValue = ["",       "",      ""     ];
        suppress = 1;
    } else if (type == "rule") {
        argList  = ["reorder"];
        defValue = ["0"      ];
        suppress = 1;
    } else if (type == "scale") {
        argList  = ["fact"];
        defValue = [""    ];
        suppress = 1;
    } else if (type == "select") {
        argList  = ["$type", "arg1", "arg2", "arg3", "arg4", "arg5", "arg6", "arg7", "arg8"];
        defValue = ["",      "",     "",     "",     "",     "",     "",     "",     ""    ];
    } else if (type == "set") {
        argList  = ["$pmtrName", "exprs"];
        defValue = ["",          ""     ];
    } else if (type == "skbeg") {
        argList  = ["x", "y", "z", "relative"];
        defValue = ["",  "",  "",  "1"       ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skcon") {
        argList  = ["$type", "index1", "index2", "$value"];
        defValue = ["",      "",       "-1",     "0"     ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skend") {
        argList  = ["wireOnly"];
        defValue = ["0"];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skvar") {
        argList  = ["$type", "valList"];
        defValue = ["",      ""       ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "solbeg") {
        argList  = ["$varList"];
        defValue = [""        ];
    } else if (type == "solcon") {
        argList  = ["$expr"];
        defValue = [""     ];
    } else if (type == "solend") {
        argList  = [];
        defValue = [];
    } else if (type == "sphere") {
        argList  = ["xcent", "ycent", "zcent", "radius"];
        defValue = ["",      "",      "",      ""      ];
        suppress = 1;
    } else if (type == "spline") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "store") {
        argList  = ["$name", "index", "keep"];
        defValue = ["",      "0" ,    "0"   ];
    } else if (type == "subtract") {
        argList  = ["$order", "index", "maxtol"];
        defValue = ["none",   "1",     "0"     ];
    } else if (type == "sweep") {
        argList  = [];
        defValue = [];
        suppress = 1;
    } else if (type == "throw") {
        argList  = ["errCode"];
        defValue = [""       ];
    } else if (type == "torus") {
        argList  = ["xcent", "ycent", "zcent", "dxaxis", "dyaxis", "dzaxis", "majorRad", "minorRad"];
        defValue = ["",      "",      "",      "",       "",       "",       "",         ""        ];
        suppress = 1;
    } else if (type == "translate") {
        argList  = ["dx", "dy", "dz"];
        defValue = ["",   "",   ""  ];
        suppress = 1;
    } else if (type == "udparg") {
        argList  = ["$primtype", "$argName1", "argValue1", "$argName2", "argValue2", "$argName3", "argValue3", "$argName4", "argValue4"];
        defValue = ["",          "",          "",          "",          "",          "",          "",          "",          ""         ];
    } else if (type == "udprim") {
        argList  = ["$primtype", "$argName1", "argValue1", "$argName2", "argValue2", "$argName3", "argValue3", "$argName4", "argValue4"];
        defValue = ["",          "",          "",          "",          "",          "",          "",          "",          ""         ];
        suppress = 1;
    } else if (type == "union") {
        argList  = ["tomark", "trimList", "maxtol"];
        defValue = ["0",      "0",        "0"     ];
    } else {
        return 1;
    }

    // set up the activity
    if (suppress) {
        if        (brch[ibrch].actv == 300) {     // Active
            editBrchForm.activity.value = 300;
        } else if (brch[ibrch].actv == 301) {     // Suppressed
            editBrchForm.activity.value = 301;
        } else {                                  // Inactive or Deferred
            editBrchForm.activity.value = 300;
        }
        editBrchForm.activity.style.display = "table-row";
    } else {
        editBrchForm.activity.style.display = "none";
    }

    // set up the number of Attributes
    if (wv.curMode != 1) {
        document.getElementById("numArgs").firstChild["data"]  = brch[ibrch].attrs.length;
    }

    // set up the arguments
    wv.numArgs = argList.length;

    if (wv.curMode != 1 && (type == "select" || type == "udparg" || type == "udprim")) {
        if        (brch[ibrch].args[0] == "" || brch[ibrch].args[0] == "$") {
            wv.numArgs = 0;
        } else if (brch[ibrch].args[1] == "" || brch[ibrch].args[1] == "$") {
            wv.numArgs = 1;
        } else if (brch[ibrch].args[2] == "" || brch[ibrch].args[2] == "$") {
            wv.numArgs = 2;
        } else if (brch[ibrch].args[3] == "" || brch[ibrch].args[3] == "$") {
            wv.numArgs = 3;
        } else if (brch[ibrch].args[4] == "" || brch[ibrch].args[4] == "$") {
            wv.numArgs = 4;
        } else if (brch[ibrch].args[5] == "" || brch[ibrch].args[5] == "$") {
            wv.numArgs = 5;
        } else if (brch[ibrch].args[6] == "" || brch[ibrch].args[6] == "$") {
            wv.numArgs = 6;
        } else if (brch[ibrch].args[7] == "" || brch[ibrch].args[7] == "$") {
            wv.numArgs = 7;
        } else if (brch[ibrch].args[8] == "" || brch[ibrch].args[8] == "$") {
            wv.numArgs = 8;
        } else if (brch[ibrch].args[0] == "" || brch[ibrch].args[0] == "$") {
            wv.numArgs = 9;
        }
    }

    if (wv.numArgs >= 1) {
        document.getElementById("argName1").parentNode.style.display = "table-row";
        document.getElementById("argName1").firstChild["data"]       = argList[0];
        if (wv.curMode == 1 || brch[ibrch].args[0] === undefined) {
            editBrchForm.argValu1.value = defValue[0];
        } else if (argList[0].charAt(0) == "$") {
            editBrchForm.argValu1.value = brch[ibrch].args[0].substr(1,brch[ibrch].args[0].length);
        } else {
            editBrchForm.argValu1.value = brch[ibrch].args[0];
        }

        wv.getFocus = editBrchForm.argValu1;
    }
    if (wv.numArgs >= 2) {
        document.getElementById("argName2").parentNode.style.display = "table-row";
        document.getElementById("argName2").firstChild["data"]       = argList[1];
        if (wv.curMode == 1 || brch[ibrch].args[1] === undefined) {
            editBrchForm.argValu2.value = defValue[1];
        } else if (argList[1].charAt(0) == "$") {
            editBrchForm.argValu2.value = brch[ibrch].args[1].substr(1,brch[ibrch].args[1].length);
        } else {
            editBrchForm.argValu2.value = brch[ibrch].args[1];
        }
    }
    if (wv.numArgs >= 3) {
        document.getElementById("argName3").parentNode.style.display = "table-row";
        document.getElementById("argName3").firstChild["data"]       = argList[2];
        if (wv.curMode == 1 || brch[ibrch].args[2] === undefined) {
            editBrchForm.argValu3.value = defValue[2];
        } else if (argList[2].charAt(0) == "$") {
            editBrchForm.argValu3.value = brch[ibrch].args[2].substr(1,brch[ibrch].args[2].length);
        } else {
            editBrchForm.argValu3.value = brch[ibrch].args[2];
        }
    }
    if (wv.numArgs >= 4) {
        document.getElementById("argName4").parentNode.style.display = "table-row";
        document.getElementById("argName4").firstChild["data"]       = argList[3];
        if (wv.curMode == 1 || brch[ibrch].args[3] === undefined) {
            editBrchForm.argValu4.value = defValue[3];
        } else if (argList[3].charAt(0) == "$") {
            editBrchForm.argValu4.value = brch[ibrch].args[3].substr(1,brch[ibrch].args[3].length);
        } else {
            editBrchForm.argValu4.value = brch[ibrch].args[3];
        }
    }
    if (wv.numArgs >= 5) {
        document.getElementById("argName5").parentNode.style.display = "table-row";
        document.getElementById("argName5").firstChild["data"]       = argList[4];
        if (wv.curMode == 1 || brch[ibrch].args[4] === undefined) {
            editBrchForm.argValu5.value = defValue[4];
        } else if (argList[4].charAt(0) == "$") {
            editBrchForm.argValu5.value = brch[ibrch].args[4].substr(1,brch[ibrch].args[4].length);
        } else {
            editBrchForm.argValu5.value = brch[ibrch].args[4];
        }
    }
    if (wv.numArgs >= 6) {
        document.getElementById("argName6").parentNode.style.display = "table-row";
        document.getElementById("argName6").firstChild["data"]       = argList[5];
        if (wv.curMode == 1 || brch[ibrch].args[5] === undefined) {
            editBrchForm.argValu6.value = defValue[5];
        } else if (argList[5].charAt(0) == "$") {
            editBrchForm.argValu6.value = brch[ibrch].args[5].substr(1,brch[ibrch].args[5].length);
        } else {
            editBrchForm.argValu6.value = brch[ibrch].args[5];
        }
    }
    if (wv.numArgs >= 7) {
        document.getElementById("argName7").parentNode.style.display = "table-row";
        document.getElementById("argName7").firstChild["data"]       = argList[6];
        if (wv.curMode == 1 || brch[ibrch].args[6] === undefined) {
            editBrchForm.argValu7.value = defValue[6];
        } else if (argList[6].charAt(0) == "$") {
            editBrchForm.argValu7.value = brch[ibrch].args[6].substr(1,brch[ibrch].args[6].length);
        } else {
            editBrchForm.argValu7.value = brch[ibrch].args[6];
        }
    }
    if (wv.numArgs >= 8) {
        document.getElementById("argName8").parentNode.style.display = "table-row";
        document.getElementById("argName8").firstChild["data"]       = argList[7];
        if (wv.curMode == 1 || brch[ibrch].args[7] === undefined) {
            editBrchForm.argValu8.value = defValue[7];
        } else if (argList[7].charAt(0) == "$") {
            editBrchForm.argValu8.value = brch[ibrch].args[7].substr(1,brch[ibrch].args[7].length);
        } else {
            editBrchForm.argValu8.value = brch[ibrch].args[7];
        }
    }
    if (wv.numArgs >= 9) {
        document.getElementById("argName9").parentNode.style.display = "table-row";
        document.getElementById("argName9").firstChild["data"]       = argList[8];
        if (wv.curMode == 1 || brch[ibrch].args[8] === undefined) {
            editBrchForm.argValu9.value = defValue[8];
        } else if (argList[8].charAt(0) == "$") {
            editBrchForm.argValu9.value = brch[ibrch].args[8].substr(1,brch[ibrch].args[8].length);
        } else {
            editBrchForm.argValu9.value = brch[ibrch].args[8];
        }
    }

    if (wv.curMode != 1) {
        if (brch[ibrch].attrs.length > 0) {
            document.getElementById("attrName1").parentNode.style.display = "table-row";
            document.getElementById("attrName1").firstChild["data"]       = brch[ibrch].attrs[0][0];
            document.getElementById("attrType1").firstChild["data"]       = brch[ibrch].attrs[0][1];
            editBrchForm.attrValu1.value                                  = brch[ibrch].attrs[0][2];
        }
        if (brch[ibrch].attrs.length > 1) {
            document.getElementById("attrName2").parentNode.style.display = "table-row";
            document.getElementById("attrName2").firstChild["data"]       = brch[ibrch].attrs[1][0];
            document.getElementById("attrType2").firstChild["data"]       = brch[ibrch].attrs[1][1];
            editBrchForm.attrValu2.value                                  = brch[ibrch].attrs[1][2];
        }
        if (brch[ibrch].attrs.length > 2) {
            document.getElementById("attrName3").parentNode.style.display = "table-row";
            document.getElementById("attrName3").firstChild["data"]       = brch[ibrch].attrs[2][0];
            document.getElementById("attrType3").firstChild["data"]       = brch[ibrch].attrs[2][1];
            editBrchForm.attrValu3.value                                  = brch[ibrch].attrs[2][2];
        }
        if (brch[ibrch].attrs.length > 3) {
            document.getElementById("attrName4").parentNode.style.display = "table-row";
            document.getElementById("attrName4").firstChild["data"]       = brch[ibrch].attrs[3][0];
            document.getElementById("attrType4").firstChild["data"]       = brch[ibrch].attrs[3][1];
            editBrchForm.attrValu4.value                                  = brch[ibrch].attrs[3][2];
        }
        if (brch[ibrch].attrs.length > 4) {
            document.getElementById("attrName5").parentNode.style.display = "table-row";
            document.getElementById("attrName5").firstChild["data"]       = brch[ibrch].attrs[4][0];
            document.getElementById("attrType5").firstChild["data"]       = brch[ibrch].attrs[4][1];
            editBrchForm.attrValu5.value                                  = brch[ibrch].attrs[4][2];
        }
        if (brch[ibrch].attrs.length > 5) {
            document.getElementById("attrName6").parentNode.style.display = "table-row";
            document.getElementById("attrName6").firstChild["data"]       = brch[ibrch].attrs[5][0];
            document.getElementById("attrType6").firstChild["data"]       = brch[ibrch].attrs[5][1];
            editBrchForm.attrValu6.value                                  = brch[ibrch].attrs[5][2];
        }
        if (brch[ibrch].attrs.length > 6) {
            document.getElementById("attrName7").parentNode.style.display = "table-row";
            document.getElementById("attrName7").firstChild["data"]       = brch[ibrch].attrs[6][0];
            document.getElementById("attrType7").firstChild["data"]       = brch[ibrch].attrs[6][1];
            editBrchForm.attrValu7.value                                  = brch[ibrch].attrs[6][2];
        }
        if (brch[ibrch].attrs.length > 7) {
            document.getElementById("attrName8").parentNode.style.display = "table-row";
            document.getElementById("attrName8").firstChild["data"]       = brch[ibrch].attrs[7][0];
            document.getElementById("attrType8").firstChild["data"]       = brch[ibrch].attrs[7][1];
            editBrchForm.attrValu8.value                                  = brch[ibrch].attrs[7][2];
        }
        if (brch[ibrch].attrs.length > 8) {
            document.getElementById("attrName9").parentNode.style.display = "table-row";
            document.getElementById("attrName9").firstChild["data"]       = brch[ibrch].attrs[8][0];
            document.getElementById("attrType9").firstChild["data"]       = brch[ibrch].attrs[8][1];
            editBrchForm.attrValu9.value                                  = brch[ibrch].attrs[8][2];
        }
    }

    return 0;
}


//
// load info into editPmtrForm
//
function setupEditPmtrForm() {
    // alert("setupEditPmtrForm()");

    var ipmtr = wv.curPmtr;
    var name  = pmtr[ipmtr].name;
    var nrow  = pmtr[ipmtr].nrow;
    var ncol  = pmtr[ipmtr].ncol;

    var editPmtrForm = document.getElementById("editPmtrForm");

    // fill in the Parameter name
    document.getElementById("pmtrName").firstChild["data"] = name;

    var pmtrTable = document.getElementById("editPmtrTable");

    // remove old table entries
    if (pmtrTable) {
        var child1 = pmtrTable.lastChild;
        while (child1) {
            var child2 = child1.lastChild;
            while (child2) {

                var child3 = child2.lastChild;
                while (child3) {
                    child2.removeChild(child3);
                    child3 = child2.lastChild;
                }
                child1.removeChild(child2);
                child2 = child1.lastChild;
            }
            pmtrTable.removeChild(child1);
            child1 = pmtrTable.lastChild;
        }
    }

    // build the table that will contain values
    for (var irow = 0; irow <= nrow; irow++) {
        var newTR = document.createElement("TR");
        pmtrTable.appendChild(newTR);

        // fill the row
        if (irow == 0) {
            var newTD = document.createElement("TD");
            newTR.appendChild(newTD);

            var newText = document.createTextNode("");
            newTD.appendChild(newText);

            for (var icol = 1; icol <= ncol; icol++) {
                newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                newText = document.createTextNode("col\u00a0"+icol);
                newTD.appendChild(newText);
            }
        } else{
            var newTD = document.createElement("TD");
            newTR.appendChild(newTD);

            var newText = document.createTextNode("row\u00a0"+irow);
            newTD.appendChild(newText);

            for (var icol = 1; icol <= ncol; icol++) {
                var indx = icol-1 + (irow-1)*pmtr[ipmtr].ncol;

                newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                var newInput = document.createElement("input");
                newInput.type  = "text";
                newInput.name  = "row"+irow+"col"+icol;
                newInput.size  = 12;
                newInput.value = pmtr[ipmtr].value[indx];
                newTD.appendChild(newInput);

                if (irow == 1 && icol == 1) {
                    wv.getFocus = newInput;
                }
            }
        }
    }

    return 0;
}


//
// send a message to the server
//
function browserToServer(text) {
    // alert("browserToServer(text="+text+")");

    if (wv.debugUI) {
        var date = new Date;
        console.log("("+date.toTimeString().substring(0,8)+") browser-->server: "+text.substring(0,40));
    }

    wv.socketUt.send(text);
}


//
// count number of Parameter value changes
//
function numberOfPmtrChanges() {
    // alert("in numberOfPmtrChanges()");

    var nchange = 0;

    var editPmtrForm = document.getElementById("editPmtrForm");

    var ipmtr = wv.curPmtr;

    // examine each of the values
    var index   = -1;
    for (var irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
        for (var icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
            index++;

            // get the new value
            var myInput = editPmtrForm["row"+irow+"col"+icol];
            var value   = myInput.value.replace(/\s/g, "");

            if (value != pmtr[ipmtr].value[index]) {
                nchange++;
            }
        }
    }

    // return the number of changes
    return nchange;
}


//
// count number of Branch changes
//
function numberOfBrchChanges() {
    // alert("in numberOfBrchChanges()");

    var nchange = 0;

    var editBrchForm = document.getElementById("editBrchForm");

    var ibrch = wv.curBrch;

    // check the name
    if (brch[ibrch].name != editBrchForm.brchName.value.replace(/\s/g, "")) {
        nchange++;
    }

    // check the activity
    if (brch[ibrch].actv != editBrchForm.activity.value) {
        nchange++;
    }

    // check the arguments
    for (var iarg = 0; iarg < wv.numArgs; iarg++) {
        if        (iarg == 0) {
            var name  = document.getElementById("argName1").firstChild["data"];
            var value = editBrchForm.            argValu1.value.replace(/\s/g, "");
        } else if (iarg == 1) {
            var name  = document.getElementById("argName2").firstChild["data"];
            var value = editBrchForm.            argValu2.value.replace(/\s/g, "");
        } else if (iarg == 2) {
            var name  = document.getElementById("argName3").firstChild["data"];
            var value = editBrchForm.            argValu3.value.replace(/\s/g, "");
        } else if (iarg == 3) {
            var name  = document.getElementById("argName4").firstChild["data"];
            var value = editBrchForm.            argValu4.value.replace(/\s/g, "");
        } else if (iarg == 4) {
            var name  = document.getElementById("argName5").firstChild["data"];
            var value = editBrchForm.            argValu5.value.replace(/\s/g, "");
        } else if (iarg == 5) {
            var name  = document.getElementById("argName6").firstChild["data"];
            var value = editBrchForm.            argValu6.value.replace(/\s/g, "");
        } else if (iarg == 6) {
            var name  = document.getElementById("argName7").firstChild["data"];
            var value = editBrchForm.            argValu7.value.replace(/\s/g, "");
        } else if (iarg == 7) {
            var name  = document.getElementById("argName8").firstChild["data"];
            var value = editBrchForm.            argValu8.value.replace(/\s/g, "");
        } else if (iarg == 8) {
            var name  = document.getElementById("argName9").firstChild["data"];
            var value = editBrchForm.            argValu9.value.replace(/\s/g, "");
        }

        var output;
        if (name.charAt(0) != "$" || value.length <= 0) {
            output =       value;
        } else {
            output = "$" + value;
        }

        if (output != brch[ibrch].args[iarg]) {
            nchange++;
        }
    }

    // check the attributes
    for (var iattr = 0; iattr < brch[ibrch].attrs.length; iattr++) {
        if        (iattr == 0) {
            var name  = document.getElementById("attrName1").firstChild["data"];
            var value = editBrchForm.            attrValu1.value.replace(/\s/g, "");
        } else if (iattr == 1) {
            var name  = document.getElementById("attrName2").firstChild["data"];
            var value = editBrchForm.            attrValu2.value.replace(/\s/g, "");
        } else if (iattr == 2) {
            var name  = document.getElementById("attrName3").firstChild["data"];
            var value = editBrchForm.            attrValu3.value.replace(/\s/g, "");
        } else if (iattr == 3) {
            var name  = document.getElementById("attrName4").firstChild["data"];
            var value = editBrchForm.            attrValu4.value.replace(/\s/g, "");
        } else if (iattr == 4) {
            var name  = document.getElementById("attrName5").firstChild["data"];
            var value = editBrchForm.            attrValu5.value.replace(/\s/g, "");
        } else if (iattr == 5) {
            var name  = document.getElementById("attrName6").firstChild["data"];
            var value = editBrchForm.            attrValu6.value.replace(/\s/g, "");
        } else if (iattr == 6) {
            var name  = document.getElementById("attrName7").firstChild["data"];
            var value = editBrchForm.            attrValu7.value.replace(/\s/g, "");
        } else if (iattr == 7) {
            var name  = document.getElementById("attrName8").firstChild["data"];
            var value = editBrchForm.            attrValu8.value.replace(/\s/g, "");
        } else if (iattr == 8) {
            var name  = document.getElementById("attrName9").firstChild["data"];
            var value = editBrchForm.            attrValu9.value.replace(/\s/g, "");
        }

        if (value != brch[ibrch].attrs[iattr][2]) {
            nchange++;
        }
    }

    // return the number of changes
    return nchange;
}


//
// unhighlight background colors in column 1 of myTree
//
function unhighlightColumn1() {
    // alert("in unhighlightColumn1()");

    for (var jnode = 0; jnode < myTree.name.length; jnode++) {
        var myElem = document.getElementById("node"+jnode+"col1");
        if (myElem) {
            if (myElem.className == "currentTD" || myElem.className == "parentTD" || myElem.className == "childTD") {
                myElem.className = undefined;
            }
        }
    }
}


//
// print an object and its contents
//
function printObject(obj) {
    var out = '';

    for (var p in obj) {
        out += p + ': ' + obj[p] + '\n';
    }

    alert(out);
}


//
// simple text formatter patterned after C's sprintf
//
function sprintf()
{

    // if no arguments, return an empty string
    var narg = arguments.length;
    if (narg == 0) {
        return "";
    }

    // otherwise, build output from input
    var format = arguments[0];
    var answer = "";
    var iarg   = 1;

    while (1) {
        var indx = format.indexOf("%");
        if (indx >= 0) {
            answer += format.substring(0, indx);
            if (iarg < narg) {
                answer += arguments[iarg++];
            } else {
                answer += "%";
            }
            format  = format.substring(indx+1);
        } else {
            answer += format;
            break;
        }
    }

    return answer;
}
