#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import urllib
import urllib.request
import requests
from lxml import etree
import random
import time
import MySid
import math
class Spider(object):
############	搜索结构
	def change(self):
		USER_AGENTS = [
			"Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; AcooBrowser; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
			"Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.0; Acoo Browser; SLCC1; .NET CLR 2.0.50727; Media Center PC 5.0; .NET CLR 3.0.04506)",
			"Mozilla/4.0 (compatible; MSIE 7.0; AOL 9.5; AOLBuild 4337.35; Windows NT 5.1; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
			"Mozilla/5.0 (Windows; U; MSIE 9.0; Windows NT 9.0; en-US)",
			"Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; Win64; x64; Trident/5.0; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 2.0.50727; Media Center PC 6.0)",
			"Mozilla/5.0 (compatible; MSIE 8.0; Windows NT 6.0; Trident/4.0; WOW64; Trident/4.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 1.0.3705; .NET CLR 1.1.4322)",
			"Mozilla/4.0 (compatible; MSIE 7.0b; Windows NT 5.2; .NET CLR 1.1.4322; .NET CLR 2.0.50727; InfoPath.2; .NET CLR 3.0.04506.30)",
			"Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN) AppleWebKit/523.15 (KHTML, like Gecko, Safari/419.3) Arora/0.3 (Change: 287 c9dfb30)",
			"Mozilla/5.0 (X11; U; Linux; en-US) AppleWebKit/527+ (KHTML, like Gecko, Safari/419.3) Arora/0.6",
			"Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.2pre) Gecko/20070215 K-Ninja/2.1.1",
			"Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9) Gecko/20080705 Firefox/3.0 Kapiko/3.0",
			"Mozilla/5.0 (X11; Linux i686; U;) Gecko/20070322 Kazehakase/0.4.5",
			"Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.8) Gecko Fedora/1.9.0.8-1.fc10 Kazehakase/0.5.6",
			"Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11",
			"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_3) AppleWebKit/535.20 (KHTML, like Gecko) Chrome/19.0.1036.7 Safari/535.20",
			"Opera/9.80 (Macintosh; Intel Mac OS X 10.6.8; U; fr) Presto/2.9.168 Version/11.52"
		]
	###############################################
		header=random.sample(USER_AGENTS,1)
		return ''.join(header)

	def __init__(self,sid,magazine,y1,y2,author,address,y3):
		self.params={
#			'product':'WOS',
#			'parentProduct':'WOS',
#			'search_mode':'GeneralSearch',
#			'SID':sid,
			'page':'1',
			'action':'changePageSize',
			'pageSize':'50'
#			'qid':qid,			
		}

		# network -> All -> AutoSave_WOS_General... -> Headers
		self.headers = {
			'Origin': 'https://apps.webofknowledge.com',
			'Referer': 'https://apps.webofknowledge.com/UA_GeneralSearch_input.do?product=UA&search_mode=GeneralSearch&SID=R1ZsJrXOFAcTqsL6uqh&preferencesSaved=',
			'User-Agent': self.change(),
			'Content-Type': 'application/x-www-form-urlencoded'
		}
		# network -> All -> AutoSave_WOS_General... -> Headers
		# magazine
		self.form_data = {
			'fieldCount': 1,
			'action': 'search',
			'product': 'WOS',
			'search_mode': 'GeneralSearch',
			'SID': sid,
			'max_field_count': 25,
			'formUpdated': 'true',
			'value(input1)': magazine,
			'value(select1)': 'SO',
			'value(hidInput1)': '',
			'limitStatus': 'collapsed',
			'ss_lemmatization': 'On',
			'ss_spellchecking': 'Suggest',
			'SinceLastVisit_UTC': '',
			'SinceLastVisit_DATE': '',
			'range': 'CUSTOM',
			'period': 'Year Range',
			'startYear': y1,
			'endYear': y2,
			'update_back2search_link_param': 'yes',
			'ssStatus': 'display:none',
			'ss_showsuggestions': 'ON',
			'ss_query_language': '',
			'ss_numDefaultGeneralSearchFields': 1,
			'rs_sort_by': 'PY.D;LD.D;SO.A;VL.D;PG.A;AU.A'
		}
		
		# author & address
		self.form_data1 = {
			'fieldCount': 2,
			'action': 'search',
			'product': 'WOS',
			'search_mode': 'GeneralSearch',
			'SID': sid,
			'max_field_count': 25,
			'formUpdated': 'true',
			'value(input1)': author,
			'value(select1)': 'AU',
			'value(hidInput1)': '',
			'value(bool_1_2)':'AND',
			'value(input2)':address,
			'value(select2)':'AD',
			'value(hidInput2)': '',
			'limitStatus': 'collapsed',
			'ss_lemmatization': 'On',
			'ss_spellchecking': 'Suggest',
			'SinceLastVisit_UTC': '',
			'SinceLastVisit_DATE': '',
			'range': 'ALL',
			'period': 'Range Selection',
			'startYear': '1986',
			'endYear': y3,
			'update_back2search_link_param': 'yes',
			'ssStatus': 'display:none',
			'ss_showsuggestions': 'ON',
			'ss_query_language': '',
			'ss_numDefaultGeneralSearchFields': 1,
			'rs_sort_by': 'PY.D;LD.D;SO.A;VL.D;PG.A;AU.A'
		}
		

		#	搜索历史??
		self.form_data2 = {
			'product': 'WOS',
			'prev_search_mode': 'CombineSearches',
			'search_mode': 'CombineSearches',
			'SID': sid,
			'action': 'remove',
			'goToPageLoc': 'SearchHistoryTableBanner',
			'currUrl': 'https://apps.webofknowledge.com/WOS_CombineSearches_input.do?SID=' + sid + '&product=WOS&search_mode=CombineSearches',
			'x': 48,
			'y': 9,
			'dSet': 1
		}
############	搜索界面，并删除
	def Search(self,url,name,s_num):
		timeout=15    
		while True:
			try:
				s1 = requests.Session()
				t_r=math.ceil(random.uniform(1,30)) 
				time.sleep(45+t_r)
				r1 = s1.post(url, data=self.form_data, headers=self.headers,timeout=timeout)
				r1.encoding = r1.apparent_encoding
				tree = etree.HTML(r1.text)
				title=tree.xpath("//value[@lang_id='']/text()")
				s_num=s_num+1
				if timeout > 60:
					print('Spider-Search:\t\n Search Page can\'t connection')
					break
				if len(title)==0:
					if s_num % 100 ==0:
						sid=re.findall(r'SID=\w+&', url)[0].replace('SID=', '').replace('&', '')
						ss=MySid.MySid(num=s_num)
						ff=ss.Sid()
						sid1=ff[0]
						s_num=ff[1]
						url=url.replace('SID=%s&'%sid,'SID=%s&'%sid1)
						sid=sid1
					timeout=timeout+5
					continue
				break
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-Search:\t\n Search Page can\'t connection')
					break
		author=re.findall('>By: <\/span><a title.*?>(.*?)[;<]',r1.text)
		url1=tree.xpath("//a[@class='smallV110 snowplow-full-record']/@href")
		url=[u.replace('/full_record.do','http://apps.webofknowledge.com/full_record.do') for u in url1]
		################ 删除不合author
		nn=[]
		for i in range(0,len(author)):	
			a0=author[i].split(',')
			a1=a0[0].strip()
			a2=a0[len(a0)-1].strip()
			if (a1.lower() in name) or (a2.lower() in name):
				nn.append(i)
		au=[]
		ti=[]
		ur=[]
		if len(nn) > 0:
			au=[author[i] for i in nn]
			ti=[title[i] for i in nn]
			ur=[url[i] for i in nn]
		s1.keep_alive = False
		return (au, ti, ur)
#################搜索合要求的页面;Elements -> 点击需要的信息###########################################	
	def Page(self,url,author,s_num):
		timeout=15
		while True:
			try:
				s1 = requests.Session()
				t_r=math.ceil(random.uniform(1,30)) 
				time.sleep(45+t_r)
				cp_page = s1.post(url, data=self.form_data, headers=self.headers,timeout=timeout)
				cp_page.encoding = cp_page.apparent_encoding
				cp_tree = etree.HTML(cp_page.text)
				title=cp_tree.xpath("//div[@class='title']/text()")
				s_num=s_num+1
				if timeout > 60:
					print('Spider-Page:\t\n Page can\'t connection')
					break
				if len(title)==0:
					if s_num % 100 ==0:
						sid=re.findall(r'SID=\w+&', url)[0].replace('SID=', '').replace('&', '')
						ss=MySid.MySid(num=s_num)
						ff=ss.Sid()
						sid1=ff[0]
						s_num=ff[1]
						url=url.replace('SID=%s&'%sid,'SID=%s&'%sid1)
						sid=sid1
					timeout=timeout+5
					continue
				break
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-Page:\t\n Page can\'t connection')
					break
		doi=cp_tree.xpath("//span[@name='doi']/text()")
		t_time=re.findall('Published:</span>[\n]?<value>(.*?)<',cp_page.text)
		paper_class=re.findall('Document Type:</span>(.*?)<',cp_page.text)
		nn1=re.findall('\(%s\)<sup><b>\[\n.*?<a href=(.*?)\]'%author,cp_page.text)
		ad=[]
		address=''
		if len(nn1)>0:
			nn2=nn1[0]
			nn3=re.findall('\(\'(.*?)\'',nn2)
			for i in nn3:
				t1=cp_tree.xpath("//a[@id='%s']/text()"%i)
				t2=t1[0]
				t3=re.compile(r'\[.*\]')
				t4=re.sub(t3,'',t2).strip()
				ad.append(t4)
			if len(ad)>1:
				address='||'.join(ad)
			else:
				address=''.join(ad)
		else :
			address=''
		cp_url1=re.findall('By:<\/span>.*?href=\"(.*?)\".*?%s'%author,cp_page.text)
		cp_url=[u.replace('/DaisyOne','http://apps.webofknowledge.com/DaisyOne') for u in cp_url1]
		ut=re.findall('dais_id=(\w+?)&',cp_url[0])
		au1=re.findall('By:<\/span><a.*?>(.*?)<\/a>',cp_page.text)
		co_author=re.findall('<\/span>(.*?)\s*\(reprint author\)',cp_page.text)
		co=sorted(set(co_author),key=co_author.index) 
		if au1[0] in co:
			e_rr=co.index(au1[0])
			email0=cp_tree.xpath("//a[@class='snowplow-author-email-addresses']/text()")
			email=email0[e_rr]
			flag=[1]
		else:
			email=""
			flag=[0]
		research=re.findall('Research Areas:</span>(.*?)<',cp_page.text)
		s1.keep_alive = False
		if len(ut)==0:
			ut=['0']		
		return (doi,t_time,paper_class,[address],[email],research,ut,flag,au1)

###############################搜索界面，并依次获取##################################################
	def Email(self,s_num,author_id,author_ab):
		root_url = 'http://apps.webofknowledge.com/WOS_GeneralSearch.do'
		timeout=15
		while True:
			try:
				t_r=math.ceil(random.uniform(1,30)) 
				time.sleep(60+t_r)
				s = requests.Session()
				r = s.post(root_url, data=self.form_data1, headers=self.headers,timeout=timeout)
				r.encoding = r.apparent_encoding
				tree = etree.HTML(r.text)
				title=tree.xpath("//value[@lang_id='']/text()")	
				s_num=s_num+1					
				if timeout > 60:
					print('Spider-Email:\t\n Email can\'t connection')
					return None			
				if len(title)==0:
					timeout=timeout+5
					continue
				break
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-Email:\t\n Email can\'t connection')
					return None
		page=tree.xpath("//a[@class='paginationNext snowplow-navigation-nextpage-top']/@href")
		p_num=tree.xpath("//span[@id='pageCount.top']/text()")			
		page_num=int(p_num[0])
		url1=tree.xpath("//a[@class='smallV110 snowplow-full-record']/@href")
		url=[u.replace('/full_record.do','http://apps.webofknowledge.com/full_record.do') for u in url1]
		if page_num > 1 :
			for i in range(2,page_num+1):
				
				str_p=page[0]
				timeout=15
				while True:
					try:
						t_r=math.ceil(random.uniform(1,30)) 
						time.sleep(45+t_r)
						ss = requests.get(str_p,timeout=timeout)
						t = etree.HTML(ss.text)          		
						title_i=t.xpath("//value[@lang_id='']/text()")
						s_num=s_num+1
						if timeout > 60:
							print('Spider-Email:\t\n Email Search Page can\'t connection')
							return None			
						if len(title_i)==0:
							if s_num % 100 ==0:
								sid=re.findall(r'SID=\w+&', str_p)[0].replace('SID=', '').replace('&', '')
								ss=MySid.MySid(num=s_num)
								ff=ss.Sid()
								sid1=ff[0]
								s_num=ff[1]
								str_p=str_p.replace('SID=%s&'%sid,'SID=%s&'%sid1)
								sid=sid1
							timeout=timeout+5
							continue
						break
					except Exception as e:
						timeout=timeout+5
						if timeout > 60:
							print('Spider-Email:\t\n Email Search Page can\'t connection')
							return None
				page=t.xpath("//a[@class='paginationNext snowplow-navigation-nextpage-top']/@href")				
				url1=t.xpath("//a[@class='smallV110 snowplow-full-record']/@href")
				url_i=[u.replace('/full_record.do','http://apps.webofknowledge.com/full_record.do') for u in url1]
				title=title+title_i
				url=url+url_i
		s.keep_alive = False
		sid=re.findall(r'SID=\w+&', url[0])[0].replace('SID=', '').replace('&', '')
		for i in range(0,len(url)):
			
			p_u=url[i]
			if s_num % 100 ==0:
				sid=re.findall(r'SID=\w+&', p_u)[0].replace('SID=', '').replace('&', '')
				ss=MySid.MySid(num=s_num)
				ff=ss.Sid()
				sid1=ff[0]
				s_num=ff[1]
				p_u=p_u.replace('SID=%s&'%sid,'SID=%s&'%sid1)
				sid=sid1
			ff=self.E_Page(url=p_u,author_ab=author_ab,s_num=s_num,sid=sid,author_id=author_id)
			s_num=ff[2]
			flag=ff[0]
			if int(flag[0])==1:
				ti=title[i]
				ee=ff[1]
				return (ti,ee,sid,author_id)
		ti=''
		ee=''
		return (ti,ee,sid,'')

#################	邮箱页面###########################################	
	def E_Page(self,url,author_ab,s_num,sid,author_id):
		timeout=15
		while True:
			try:
				s1 = requests.Session()
				t_r=math.ceil(random.uniform(1,30)) 
				time.sleep(45+t_r)
				cp_page = s1.post(url, data=self.form_data, headers=self.headers,timeout=timeout)
				cp_page.encoding = cp_page.apparent_encoding
				cp_tree = etree.HTML(cp_page.text)
				title=cp_tree.xpath("//div[@class='title']/text()")
				s_num=s_num+1
				if timeout > 60:
					print('Spider-E_Page:\t\n E_Page can\'t connection\t\n'+'URL:'+url)
					break
				if len(title)==0:
					timeout=timeout+5
					continue
				break
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-E_Page:\t\n E_Page can\'t connection')
					break
		co_author=re.findall('<\/span>(.*?)\s*\(reprint author\)',cp_page.text)
		co=sorted(set(co_author),key=co_author.index)
		flag=0
		email=''
		if author_ab in co:
			if author_id=='0':
				e_rr=co.index(author_ab)
				email0=cp_tree.xpath("//a[@class='snowplow-author-email-addresses']/text()")
				email=email0[e_rr]

				flag=1				
			if author_id != '0':
				hr=cp_tree.xpath("//a[@title='Find more records by this author']/@href")
				hr1='|'.join(hr)			
				hr2=re.findall('&ut=(%s)'%author_id,hr1)
				if len(hr2)>0:
					e_rr=co.index(author_ab)
					email0=cp_tree.xpath("//a[@class='snowplow-author-email-addresses']/text()")
					email=email0[e_rr]

					flag=1
		else:
			email=""
			flag=0
		s1.keep_alive = False
		return ([flag],email,s_num,sid)
######################## 删除历史############################################
	def Delete_History(self,Spider):
		print('Delete History')
		murl='http://apps.webofknowledge.com/WOS_CombineSearches.do'
		s = requests.Session()
		timeout=15
		while True:
			try:
				s.post(murl, data=Spider.form_data2, headers=Spider.headers,timeout=timeout)
				if timeout > 60:
					print('Spider-Delete History:\t\n E_Page can\'t connection')
					break
				print('Successed!')
				return 'Delete History'		
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-Delete History:\t\n E_Page can\'t connection')
					break

	def Max_Page(self,s_num):
		root_url = 'http://apps.webofknowledge.com/WOS_GeneralSearch.do'
		timeout=15
		while True:
			try:
				s1 = requests.Session()
				t_r=math.ceil(random.uniform(1,30)) 
				time.sleep(45+t_r)
				r1 = s1.post(root_url, data=self.form_data, headers=self.headers,timeout=timeout)
				r1.encoding = r1.apparent_encoding
				tree = etree.HTML(r1.text)
				p_num=tree.xpath("//span[@id='pageCount.top']/text()")
				s_num=s_num+1
				if timeout > 60:
					print('Spider-Max:\t\n Max Search Page can\'t connection')
					break
				if len(p_num)==0:
					timeout=timeout+5
					continue
				break
			except Exception as e:
				timeout=timeout+5
				if timeout > 60:
					print('Spider-Max:\t\n Max Search Page can\'t connection')
					break
		page=tree.xpath("//a[@class='paginationNext snowplow-navigation-nextpage-top']/@href")
		reobj=re.compile(":\d{3}")
		list_p=reobj.sub("",page[0])
		str_p=list_p.replace('&&update_back2search_link_param=yes&page=2','').replace('parentQid=&','')	
		url=str_p+'&'+urllib.parse.urlencode(self.params) ## 50页第一次
		r2 = s1.post(url)
		return (p_num,url,list_p)


