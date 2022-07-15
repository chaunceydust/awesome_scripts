#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import requests
import time
from multiprocessing import Process

class MySid(Process):
	def __init__(self,num):
		self.num=num

	def Sid(self):
		root = 'http://www.webofknowledge.com/'
		timeout=15
		while True:
			try:
				time.sleep(45)
				s = requests.get(root,timeout=timeout)
				self.num=self.num+1
				if self.num % 100 == 0:
					self.num=1
				break
			except Exception as e:
				print("Sid error :\n\t"+str(e))
				timeout=timeout+5
				if timeout > 60:
					print('Sid-MySid:\t\n can\'t connection')
					break
		sid = re.findall(r'SID=\w+&', s.url)[0].replace('SID=', '').replace('&', '')
		return (sid,self.num)
		


