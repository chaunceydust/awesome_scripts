#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import MyThread
import MySid
import sys
import pandas as pd
import MyMax
################	main	##########################

if __name__ == "__main__":
	################	input	####################
	y1=int(sys.argv[1])
	y2=int(sys.argv[2])
	magazine_input=str(sys.argv[3])
	name_input=str(sys.argv[4])
	page_files=str(sys.argv[5])
	star_i=int(str(sys.argv[6]))
	save=str(sys.argv[7])
	################################################
	tt=save+'/test_Page.txt'
	with open(tt,"w") as f :
		f.write("这是个测试！")  #这句话自带文件关闭功能，不需要再写f.close()

	
	name=[]	
	for n in open(name_input):
		n1=n.strip('\n')
		name.append(n1)
	########	按杂志搜索
	s_num=0
	ss=MySid.MySid(num=s_num)
	ff=ss.Sid()
	sid=ff[0]
	s_num=ff[1]
	csvfile = open(page_files,'r')									
	title= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Title"].values
	author= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Author"].values
	mm= pd.read_csv(filepath_or_buffer = page_files ,sep = ',')["Magazine"].values
	page_url= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["URL"].values
	csvfile.close()
	#######	设置最大50页
	myMax=MyMax.MyMax(sid=sid,y1=y1,y2=y2,y3='',author=name,address='',s_num=s_num,page_url=page_url)
	page_url=myMax.Max(magazine_input,mm)

##################### 	提取网页内容	#########################
	file_result=save+'/result_Page.csv'
	out1 = open(file_result,'w',newline='')
	result1=csv.writer(out1,dialect='excel')
	result1.writerow(['Title','DOI','Time','Class','Magazine','Author','Author_Ab','Address','Email','Research Direction',"Author ID",'Flag'])
	out1.close()	
	out1 = open(file_result,'a',newline='')
	result1=csv.writer(out1,dialect='excel')

	for i in range(star_i-1,len(page_url)):
		print("Page:"+str(i+1)+'/'+str(len(page_url)))
		p_u=page_url[i]
		########### 1小时换1次sid
		if s_num % 100 ==0:
			ss=MySid.MySid(num=s_num)
			ff=ss.Sid()
			sid1=ff[0]
			s_num=ff[1]
			p_u=p_u.replace('SID=%s&'%sid,'SID=%s&'%sid1)
			myMax=MyMax.MyMax(sid=sid1,y1=y1,y2=y2,y3='',author=name,address='',s_num=s_num,page_url=page_url)
			page_url=myMax.Max(magazine_input,mm)
			sid=sid1
		#########################################
		tt = MyThread.MyThread(sid=sid,magazine=mm[i],y1=y1,y2=y2,author=author[i],address='',y3='',url=p_u,s_num=s_num,title=title[i])
		ffp=tt.Page()
		s_num=s_num+1
		if len(ffp[6])==0:
			result1.writerow([tt.title,ffp[0][0],ffp[1][0],ffp[2][0],tt.magazine,tt.author,ffp[8][0],ffp[3][0],ffp[4][0],ffp[5][0],'',ffp[7][0]])	
		else:
			result1.writerow([tt.title,ffp[0][0],ffp[1][0],ffp[2][0],tt.magazine,tt.author,ffp[8][0],ffp[3][0],ffp[4][0],ffp[5][0],ffp[6][0],ffp[7][0]])	
	out1.close()		
