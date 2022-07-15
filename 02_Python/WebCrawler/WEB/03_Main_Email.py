#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import csv
import Spider
import MyThread
import MySid
import sys
import pandas as pd
################	main	##########################

if __name__ == "__main__":
	################	input	####################
	page_files=str(sys.argv[1])
	star_i=int(str(sys.argv[2]))
	save=str(sys.argv[3])
	################################################
	tt=save+'/test_Email.txt'
	with open(tt,"w") as f :
		f.write("这是个测试！")  #这句话自带文件关闭功能，不需要再写f.close()

	today=datetime.datetime.now()
	y3=today.year
	########	按杂志搜索
	csvfile = open(page_files,'r')									
	title= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Title"].values
	author= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Author"].values
	mm= pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Magazine"].values
	doi=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["DOI"].values
	p_time=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Time"].values
	paper_class=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Class"].values
	address=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Address"].values
	email=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Email"].values
	research=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Research Direction"].values
	author_id=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Author ID"].values
	flag_email=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Flag"].values
	author_ab=pd.read_csv(filepath_or_buffer = page_files, sep = ',')["Author_Ab"].values
	csvfile.close()
	reference_ti=[]
	au_id2=[]
	for i in range(0,len(title)):
		reference_ti=reference_ti+['']
		au_id2=au_id2+['']
		author_id[i-1]=int(author_id[i-1])
	s_num=0
#######################	搜索作者与地址，获取邮箱
	file_result=save+'/result_Email.csv'
	out = open(file_result,'w',newline='')
	result=csv.writer(out,dialect='excel')
	result.writerow(['Title','DOI','Time','Class','Magazine','Author','Author_Ab','Address','Email','Research Direction',"Author ID","Reference title",'Author ID2'])
	out.close()	
	out = open(file_result,'a',newline='')
	result=csv.writer(out,dialect='excel')
	ii_num=0	#记录搜索次数
	for i in range(star_i-1,len(flag_email)):
		if flag_email[i]==0:
			print("********************************\t\nEmail:"+str(i+1)+'/'+str(len(flag_email))+' ,Author:'+author[i]+' ,Search num:'+str(ii_num))
			if type(address[i])==float:
				result.writerow([title[i],doi[i],p_time[i],paper_class[i],mm[i],author[i],author_ab[i],address[i],email[i],
							research[i],author_id[i],reference_ti[i],au_id2[i]])
				continue
			else:
				ad=address[i].split('||')
				for j in ad:
					if s_num % 100 ==0:
						ss=MySid.MySid(num=s_num)
						ff=ss.Sid()
						sid=ff[0]
						s_num=ff[1]
					if ii_num % 100==0: ######最高200,150次清空
						my1=Spider.Spider(sid=sid,magazine=mm[i],y1='',y2='',author='',address='',y3='')
						ffd=my1.Delete_History(my1)
					ii_num=ii_num+1
					tt= MyThread.MyThread(sid=sid,magazine=mm[i],y1='',y2='',author=author[i],author_ab=author_ab[i],address=address[i],y3=y3,url='',
						s_num=s_num,title=title[i],reference_ti=reference_ti[i],email=email[i],doi=doi[i],p_time=p_time[i],
						paper_class=paper_class[i],research=research[i],author_id=str(author_id[i]),i=i,address2=j)
		if flag_email[i]==1:
			result.writerow([title[i],doi[i],p_time[i],paper_class[i],mm[i],author[i],author_ab[i],address[i],email[i],
							research[i],author_id[i],reference_ti[i],au_id2[i]])
			continue
		if type(email[i])==float or len(email[i])==0:
			ffe=tt.Email()
			sid=ffe[2]
			email[i]=ffe[1]
			reference_ti[i]=ffe[0]
			au_id2[i]=ffe[3]	
		else:
			continue	
		result.writerow([title[i],doi[i],p_time[i],paper_class[i],mm[i],author[i],author_ab[i],address[i],email[i],
						research[i],author_id[i],reference_ti[i],au_id2[i]])	
	out.close()
	
###############################################

