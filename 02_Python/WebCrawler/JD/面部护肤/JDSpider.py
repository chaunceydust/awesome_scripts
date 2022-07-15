#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import requests
import random
from lxml import etree
import time
import math
import Items
from scrapy import Spider
class JDSpider(Spider):
	def __init__(self, url,key='',class0='',class_id='',category='',items_id=''):
		self.url=url
		self.key=key
		self.class0=class0
		self.class_id=class_id
		self.items_id=items_id
		self.category=category
		self.headers = {
			'Accept': '*/*',
			'Accept-Language': 'zh-CN,zh;q=0.9',
			'User-Agent': self.change(),
			'Connection': 'keep-alive',
			'Referer': 'https://www.jd.com/allSort.aspx'
		}
######################	浏览器	########################################
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
		header=random.sample(USER_AGENTS,1)
		return ''.join(header)

######################	获取分类	####################################	
	def parse_category(self):
		timeout=6   
		while True:
			try:
				s1 = requests.Session()
				
				
				r1 = s1.post(url=self.url,headers=self.headers)
				r1.encoding = r1.apparent_encoding
				tree = etree.HTML(r1.text)
				texts = tree.xpath('//dl[@class="clearfix"]')
				items=[]
				for text in texts:
					tt= etree.tostring(text,encoding="utf-8")
					ss=tt.decode('utf-8')
					items_class=re.findall(r'<a href=".*?" target="_blank">(.*?)</a>', ss)
					tmp=list(set([items_class[0]]).intersection(set(self.key))) ###求交集
					if len(tmp)==1:
						print(tmp[0])
						items_url=re.findall(r'<a href="(.*?)" target="_blank">.*?</a>', ss)
						for i in range(1,len(items_class)):
							categoriesItem=Items.CategoriesItem()
							categoriesItem['category']=items_class[i]
							categoriesItem['url']='https:'+items_url[i]
							categoriesItem['class_id']=re.findall('cat=(.*)',items_url[i])[0]
							categoriesItem['_class']=tmp[0]
							items=items+[categoriesItem]
				return (items)
				break
			except Exception as e:
				print("JDSpider-parse_category :"+str(e))
				timeout=timeout+5
				if timeout > 60:
					print('JDSpider-parse_category timeout,finish!')
					break
	

######################	获取产品	####################################
	def parse_list(self): #分别获得商品的地址和下一页地址"""
		timeout=6   
		while True:
			try:
				s1 = requests.Session()				
				r1 = s1.post(url=self.url,headers=self.headers)
				r1.encoding = r1.apparent_encoding
				tree = etree.HTML(r1.text)
				page=tree.xpath('//span[@class="p-skip"]/em/b/text()')
				if len(page)==0:
					page=['1']
				uu=[]
				cc=[]
				for i in range(0,int(page[0])):
					print('parse_list:'+str(i)+' - '+page[0])
					if i ==0:
						cc0=self.url.split('=')[1].split('&')[0]
						uu1 = tree.xpath('//div[@class="p-img"]/a/@href')
						uu2=['https:'+u for u in uu1]
						uu=uu+uu2
						if int(page[0])!=1:
							url='https://list.jd.com'+tree.xpath('//a[@class="pn-next"]/@href')[0]
						cc=cc+len(uu2)*[cc0]
					else :
						timeout=6 
						while True:
							try:
								r1 = s1.post(url=url,headers=self.headers)
								r1.encoding = r1.apparent_encoding
								tree = etree.HTML(r1.text)
								page1=tree.xpath('//span[@class="p-skip"]/em/b/text()')
								if len(page1)>0:
									break
							except Exception as e:
								print("JDSpider-parse_list :"+str(e))
								timeout=timeout+5
								if timeout > 60:
									print('JDSpider-parse_list timeout,finish!')
									break
						cc0=url.split('=')[1].split('&')[0]
						uu1 = tree.xpath('//div[@class="p-img"]/a/@href')
						uu2=['https:'+u for u in uu1]
						uu=uu+uu2
						cc=cc+len(uu2)*[cc0]
						if page1[0] != page[0]:
							url='https://list.jd.com'+tree.xpath('//a[@class="pn-next"]/@href')[0]
				print('JDSpider-parse_list : URL_Length='+str(len(uu)))
				return (uu)
				break
			except Exception as e:
				print("JDSpider-parse_list :"+str(e))
				timeout=timeout+5
				if timeout > 60:
					print('JDSpider-parse_list timeout,finish!')
					break		

	def parse_product(self):	#商品页获取title,price,product_id
		price_url = 'https://p.3.cn/prices/mgets?skuIds=J_'
		favourable_url = 'https://cd.jd.com/promotion/v2?skuId=%s&area=1_72_2799_0&shopId=%s&venderId=%s&cat=%s'
		Base_url = 'https://list.jd.com'
		comment_url = 'https://club.jd.com/comment/productPageComments.action?productId=%s&score=0&sortType=5&page=%s&pageSize=10'
		timeout=6   
		while True:
			try:
				s1 = requests.Session()
				r1 = s1.post(url=self.url,headers=self.headers)
				r1.encoding = r1.apparent_encoding
				tree = etree.HTML(r1.text)
				########### 商铺信息
				ids = re.findall(r"venderId:(.*?),\s.*?shopId:'(.*?)'", r1.text)
				if not ids:
					ids = re.findall(r"venderId:(.*?),\s.*?shopId:(.*?),", r1.text)
				vender_id = ids[0][0]
				shop_id = ids[0][1]
				shopItem = Items.ShopItem()
				shopItem['shopId'] = shop_id
				shopItem['venderId'] = vender_id
				shopItem['url1'] = 'http://mall.jd.com/index-%s.html' % (shop_id)
				try:
					shopItem['url2'] = 'https:' + tree.xpath('//ul[@class="parameter2 p-parameter-list"]/li/a/@href')[0]
				except:
					shopItem['url2'] = shopItem['url1']
				name = ''
				if shop_id == '0':
					name = '京东自营'
				else:
					try:
						name = tree.xpath('//ul[@class="parameter2 p-parameter-list"]/li/a//text()')[0]
					except:
						try:
							name = tree.xpath('//div[@class="name"]/a//text()')[0]
						except:
							try:
								name = tree.xpath('//div[@class="shopName"]/strong/span/a//text()')[0]
							except:
								try:
									name = tree.xpath('//div[@class="seller-infor"]/a//text()')[0]
								except:
									name = u'京东自营'
				shopItem['name'] = name
				shopItem['_id'] = name
		#############	抓取商品
				productsItem = Items.ProductsItem()
				productsItem['shopId'] = shop_id
				productsItem['category'] = self.category
				productsItem['_class'] = self.class0
				productsItem['class_id'] = self.class_id
				productsItem['items_id'] = self.items_id
				try:
					title = tree.xpath('//div[@class="sku-name"]/text()')[0].replace(u"\xa0", "").strip()
				except Exception as e:
					title = tree.xpath('//div[@id="name"]/h1/text()')[0]
				productsItem['name'] = title
				product_id = self.url.split('/')[-1][:-5]
				productsItem['_id'] = product_id
				productsItem['items_url'] = self.url
      		  # description
				desc = tree.xpath('//ul[@class="parameter2 p-parameter-list"]//text()')
				productsItem['description'] = ';'.join(i.strip() for i in desc)
	
     		   # price
				response = requests.get(url=price_url + product_id)
				price_json = response.json()
				productsItem['reallyPrice'] = price_json[0]['p']
				productsItem['originalPrice'] = price_json[0]['m']

        # 优惠
				res_url = favourable_url % (product_id, shop_id, vender_id, self.class_id.replace(',', '%2c'))
     	   # print(res_url)
				response = requests.get(res_url)
				fav_data = response.json()
				if fav_data['skuCoupon']:
					desc1 = []
					for item in fav_data['skuCoupon']:
						start_time = item['beginTime']
						end_time = item['endTime']
						time_dec = item['timeDesc']
						fav_price = item['quota']
						fav_count = item['discount']
						fav_time = item['addDays']
						desc1.append(u'有效期%s至%s,满%s减%s' % (start_time, end_time, fav_price, fav_count))
					productsItem['favourableDesc1'] = ';'.join(desc1)
				else :
					productsItem['favourableDesc1'] =''

				if fav_data['prom'] and fav_data['prom']['pickOneTag']:
					desc2 = []
					for item in fav_data['prom']['pickOneTag']:
						desc2.append(item['content'])
					productsItem['favourableDesc2'] = ';'.join(desc2)	
				else :
					productsItem['favourableDesc2'] =''	
        # 评论数
				cc_url=comment_url % (product_id, '0')			
				response = requests.get(cc_url)
				cc_data = response.json()
				productsItem["goodRateShow"]=cc_data['productCommentSummary']["goodRateShow"]
				productsItem["poorRateShow"]=cc_data['productCommentSummary']["poorRateShow"]
				productsItem["generalRateShow"]=cc_data['productCommentSummary']["generalRateShow"]
				productsItem["goodCount"]=cc_data['productCommentSummary']["goodCount"]
				productsItem["generalCount"]=cc_data['productCommentSummary']["generalCount"]
				productsItem["poorCount"]=cc_data['productCommentSummary']["poorCount"]
				productsItem["commentCount"]=cc_data['productCommentSummary']["commentCount"]
				productsItem["defaultGoodCount"]=cc_data['productCommentSummary']["defaultGoodCount"]	
				productsItem["averageScore"]=cc_data['productCommentSummary']["averageScore"]		

				hotCommentTagItems=[]
				for hotComment in cc_data['hotCommentTagStatistics']:
					hotCommentTagItem = Items.HotCommentTagItem()
					hotCommentTagItem['name'] = hotComment.get('name')
					hotCommentTagItem['count'] = hotComment.get('count')
					hotCommentTagItem['_class'] = self.class0
					hotCommentTagItem['category'] = self.category
					hotCommentTagItem['productId'] = product_id
					hotCommentTagItems=hotCommentTagItems+[hotCommentTagItem]
				return (shopItem,productsItem,hotCommentTagItems)
			except Exception as e:
				print("JDSpider-parse_product :"+str(e))
				timeout=timeout+5
				if timeout > 60:
					print('JDSpider-parse_product timeout,finish!')
					break
