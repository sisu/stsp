#!/usr/bin/python

from socket import *
from threading import Thread
import subprocess

class clithread(Thread):
	def __init__(self,cli):
		Thread.__init__(self)
		self.cli = cli
	def run(self):
		print "asd"
		s = ""
		while True:
			dat = self.cli.recv(1024)
			if dat==None: break
			a = dat.find(';')
			if (a>=0):
				s += dat[0:a]
				break
			s += dat

		print s
		p = subprocess.Popen(['./solve', 't3-2.gr', '-s'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#		p = subprocess.Popen(['cat'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		p.stdin.write(s)
		p.stdin.close()
		out = p.stdout.read()
		p.wait()
		print "out", out
		self.cli.send(out)
		self.cli.close()


s = socket(AF_INET, SOCK_STREAM)
s.setsockopt(SOL_SOCKET, SO_REUSEADDR, 1)
s.bind(('', 34343))
s.listen(10)

while True:
	print "loop"
	(cli,addr) = s.accept()
	ct = clithread(cli)
	ct.start()
