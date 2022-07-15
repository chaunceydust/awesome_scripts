#!/usr/bin/env python
# -*- coding: utf-8 -*-
from multiprocessing import Process
import Spider

		
class MyThread(Process):
	def __init__(self, sid,magazine,y1,y2,author,address,y3,url,author_ab='',s_num='',title='',reference_ti='',email='',doi='',p_time='',
				paper_class='',research='',author_id='',i='',address2=''):
		Process.__init__(self)
		self.sid = sid
		self.magazine = magazine
		self.y1=y1
		self.y2=y2
		self.y3=y3
		self.author=author
		self.address=address
		self.url=url
		self.author_ab=author_ab
		self.num=s_num
		self.title=title
		self.reference_ti=reference_ti
		self.email=email
		self.doi=doi
		self.p_time=p_time
		self.paper_class=paper_class
		self.research=research
		self.author_id=author_id
		self.i=i
		self.address2=address2


	def Search(self):
		ff= Spider.Spider(sid=self.sid,magazine=self.magazine,y1=self.y1,y2=self.y2,author=self.author,address=self.address,
		y3=self.y3).Search(url=self.url,name=self.author,s_num=self.num)
		return ff    

	def Page(self):
		ff= Spider.Spider(sid=self.sid,magazine=self.magazine,y1=self.y1,y2=self.y2,author=self.author,address=self.address,
		y3=self.y3).Page(url=self.url,author=self.author,s_num=self.num)
		return ff

	def Email(self):
		ff= Spider.Spider(sid=self.sid,magazine=self.magazine,y1=self.y1,y2=self.y2,author=self.author,address=self.address2,
		y3=self.y3).Email(s_num=self.num,author_id=self.author_id,author_ab=self.author_ab)
		return ff





