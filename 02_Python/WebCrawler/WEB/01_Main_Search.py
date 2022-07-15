#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import csv
import Spider
import MyThread
import MySid
import math
import sys
################	main	##########################

if __name__ == "__main__":
	################	input	####################
	y1=int(sys.argv[1])
	y2=int(sys.argv[2])
	magazine_input=str(sys.argv[3])
	name_input=str(sys.argv[4])
	save=str(sys.argv[5])
	star_i=int(str(sys.argv[6]))
	################################################
	tt=save+'/test_Search.txt'
	with open(tt,"w") as f :
		f.write("这是个测试！")  #这句话自带文件关闭功能，不需要再写f.close()


	name=[]	
	for n in open(name_input):
		n1=n.strip('\n')
		name.append(n1)

	#####	
	today=datetime.datetime.now()
	y3=today.year
	########	按杂志搜索
	s_num=0
	ss=MySid.MySid(num=s_num)
	ff=ss.Sid()
	sid=ff[0]
	s_num=ff[1]
	print('**************	Search	*********************')
	file_result=save+'/result_Search.csv'
	out0 = open(file_result,'w',newline='')
	result0=csv.writer(out0,dialect='excel')
	result0.writerow(['Title','Author','URL','Magazine'])
	out0.close()	
	out0 = open(file_result,'a',newline='')
	result0=csv.writer(out0,dialect='excel')
	#######	设置最大50页
	for m in open(magazine_input): # 按行读，参数
		m1=m.strip('\n')	#去除换行符
		print(m1)
		#############	第一次搜索，然后使每页显示50条
		my1=Spider.Spider(sid=sid,magazine=m1,y1=y1,y2=y2,author='',address='',y3='')
		ffm=my1.Max_Page(s_num)
		p_num=''.join(ffm[0])
		url=''.join(ffm[1])
		list_p=''.join(ffm[2])
		############# 多线程搜索 ###############
		page_num=int(p_num)
		max_page=math.ceil(page_num/5)
		
		#####################	提取搜索页内容	#######################
		for i in range(star_i-1,max_page):
			print(m1+"-"+str(i+1)+'/'+str(max_page))
			#################### 100次换1次sid
			if s_num % 100 ==0:
				ss=MySid.MySid(num=s_num)
				ff=ss.Sid()
				sid1=ff[0]
				s_num=ff[1]
			#########################################
			if i==0:
				if s_num % 100 ==0:
					url=url.replace('SID=%s&'%sid,'SID=%s&'%sid1)
					sid=sid1
				tt = MyThread.MyThread(sid=sid,magazine=m1,y1=y1,y2=y2,author=name,address='',y3=y3,url=url,s_num=s_num)
			else :
				if s_num % 100 ==0:
					list_p=list_p.replace('SID=%s&'%sid,'SID=%s&'%sid1)
					sid=sid1
				uu=list_p.replace('&page=2','&page=%s'%i)
				tt= MyThread.MyThread(sid=sid,magazine=m1,y1=y1,y2=y2,author=name,address='',y3=y3,url=uu,s_num=s_num)
			fft=tt.Search()
			s_num=s_num+1
			for ir in range(0,len(fft[0])):
				result0.writerow([fft[1][ir],fft[0][ir],fft[2][ir],m1])			
	out0.close()
	

