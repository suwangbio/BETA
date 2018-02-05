# print current time and information on screen
import subprocess
from subprocess import call as subpcall
import sys, os, time
import logging

error   = logging.critical		# function alias
warn    = logging.warning
def Info(infoStr):
    print "[%s] %s" %(time.strftime('%H:%M:%S'), infoStr)
    
def run_cmd(command):
    subpcall (command, shell = True)

def run_pip(command):
    command = subprocess.Popen(command,
                               shell = True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               )
