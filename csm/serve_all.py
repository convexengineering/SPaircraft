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
    os.remove("d82_000.egads")
except OSError:
    pass
Popen(["python", os.sep.join(["..","server.py"])])
wait_until_file_exists("d82-0.csm")
Popen(["serveCSM", "d82-0.csm"])
wait_until_file_exists("d82_000.egads")
sleep(5)
webbrowser.open_new(os.sep.join(["..", "ESP", "ESP-localhost7681.html"]))
