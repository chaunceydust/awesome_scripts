#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Spider
import MyThread
import re
from multiprocessing import Process
class MyMax(Process):
	def __init__(self,sid,y1,y2,y3,author,address,s_num,page_url):
		self.sid=sid
		self.y1=y1
		self.y2=y2
		self.author=author
		self.address=address
		self.y3=y3
		self.s_num=s_num
		self.page_url=page_url

	def Max(self,magazine_input,mm):
		for m in open(magazine_input): # 按行读，参数
			m1=m.strip('\n')	#去除换行符
			#############	第一次搜索，然后使每页显示50条
			my1=Spider.Spider(sid=self.sid,magazine=m1,y1=self.y1,y2=self.y2,author=self.author,address=self.address,y3=self.y3)
			ffm=my1.Max_Page(self.s_num)
			url=''.join(ffm[1])
			tt = MyThread.MyThread(sid=self.sid,magazine=m1,y1=self.y1,y2=self.y2,author=self.author,address=self.address,y3=self.y3,url=url,s_num=self.s_num)
			fft=tt.Search()
			self.s_num=self.s_num+1
			page_url=self.page_url
			for j in range(0,len(mm)):
				if mm[j]==m1:
					uu=self.page_url[j]
					tt1=re.findall(r'qid=\d+&SID=\w+&', uu)[0]
					tt2=re.findall(r'qid=\d+&SID=\w+&', url)[0]
					uu1=uu.replace('%s'%tt1,'%s'%tt2)
					page_url[j]=uu1
		return page_url

