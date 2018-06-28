import os
import sys
import webbrowser
from subprocess import Popen, call
from time import sleep


def wait_until_file_exists(filename):
    while not os.path.exists(filename):
        sleep(1)
    print "found", filename


try:
    os.remove("d82-0.csm")
    #.remove("d82_000.egads")
except OSError:
    pass
Popen(["python", "server.py"])
wait_until_file_exists("ESP/d82-0.csm")
Popen(["serveCSM","ESP/d82-0.csm"])
#wait_until_file_exists("ESP/d82_000.egads")
sleep(5)
webbrowser.open_new(os.sep.join(["..", "ESP", "ESP-localhost7681.html"]))
