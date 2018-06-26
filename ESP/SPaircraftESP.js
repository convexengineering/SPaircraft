var url = "localhost:8000" //prompt("Enter hostname:port for the gp server", "localhost:8000")

window.gp = {
  optcount: 0,

  dom: {},

  sol: {},

  esp: {build_outdated: false,
        build: function() {
          postMessage("Attempting to open \""+wv.filename+"\" ...");
          browserToServer("open|"+wv.filename+"|");
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
        },
        update: function() {
            gp.dom.buildButton.disabled = false
            window.oldactivateBuildButton()
          },
        },

  awaiting_response: false,
  last_sent: null,
  sendpmtrs: function() {
    gp.dom.optimizeButton.innerText = "Optimizing..."
    gp.dom.optimizeButton.style.backgroundColor = "#DA70D6"
    gp.optcount++

    for (var i=0; i < pmtr.length; i++) {
      console.log(i, pmtr[i].name, pmtr[i].value[0], gp.sol[pmtr[i].name])
      gp.sol[pmtr[i].name] = pmtr[i].value[0]
    }

    console.log(gp.sol)
    gp.websocket.send(JSON.stringify(gp.sol))
    // gp.websocket.send("sol")
  },

  websocket: new WebSocket("ws://"+url+"/")
}

// gp.websocket.onopen = gp.sendpmtrs
gp.websocket.onmessage = function(evt) {
  data = JSON.parse(evt.data);
  console.log("Data received:", data)
  postMessage("GP: " + data.msg)
  console.log(data.status)
  if (data.status == "optimal") {
    gp.dom.optimizeButton.innerText = "Optimized"
    gp.dom.optimizeButton.style.backgroundColor = null
    gp.esp.update()
  } else {
    gp.dom.optimizeButton.innerText = "Error"
    gp.dom.optimizeButton.style.backgroundColor = "#FF3F3F"
  }
}

window.oldactivateBuildButton = activateBuildButton
window.activateBuildButton = function() {
  gp.dom.optimizeButton.innerText = "Press to Optimize"
  gp.dom.optimizeButton.style.backgroundColor = "#3FFF3F";
}

gp.dom.buttonForm = document.getElementById("butnfrm")
gp.dom.buildButton = document.getElementById("buildButton")
gp.dom.buildButton.disabled = true
gp.dom.buildButton.onclick = gp.esp.build

if (!document.getElementById("optButton"))
  gp.dom.optimizeButton = document.createElement("button")
else
  gp.dom.optimizeButton = document.getElementById("optButton")
gp.dom.optimizeButton.id = "optButton"
gp.dom.optimizeButton.type = "button"
gp.dom.optimizeButton.innerText = "Optimized"
gp.dom.optimizeButton.onclick = gp.sendpmtrs
gp.dom.buttonForm.insertBefore(gp.dom.optimizeButton, gp.dom.buildButton)
