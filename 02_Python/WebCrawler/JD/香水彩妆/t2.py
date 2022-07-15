#!/usr/bin/env python
# -*- coding: utf-8 -*-
import JDSpider
import Items
import re
################	main	##########################

if __name__ == "__main__":
	key_word = ['book', 'e', 'channel', 'mvd', 'list']
	Base_url = 'https://list.jd.com'
	price_url = 'https://p.3.cn/prices/mgets?skuIds=J_'
	comment_url = 'https://club.jd.com/comment/productPageComments.action?productId=%s&score=0&sortType=5&page=%s&pageSize=10'
	favourable_url = 'https://cd.jd.com/promotion/v2?skuId=%s&area=1_72_2799_0&shopId=%s&venderId=%s&cat=%s'
	url='https://www.jd.com/allSort.aspx'
	key=['清洁用品','面部护肤','口腔护理','女性护理','洗发护发','香水彩妆']
#######################
	cc=[]
	files=open('Parse_Category_new.txt')
	i=0
	for line in files.readlines():
		if i==0:
			i=2
			continue
		cutline=line.strip().split('\t')
		c0=Items.CategoriesItem()
		c0['category']=cutline[0]
		c0['url']=cutline[1]
		c0['_class']=cutline[2]
		c0['class_id']=cutline[3]
		cc=cc+[c0]
		
###########		抓取产品url
	items=[]
	for i in range(0,len(cc)):
		n1=	cc[i]['category']	
		u1=cc[i]['url']
		c1=cc[i]['_class']
		i1=cc[i]['class_id']
		tt=JDSpider.JDSpider(url=u1,category=n1,class0=c1,class_id=i1)
		print('Product_list:'+str(i)+' - '+str(len(cc))+'\t\n'+u1)

		fft=tt.parse_list()
		for j in range(0,len(fft)):
			it=Items.ProductsItem()
			it['_class']=c1
			it['class_id']=i1
			it['category']=n1
			it['items_url']=fft[j]
			it['class_url']=u1
			it['items_id']=re.findall(r'(\d+)',fft[j])[0]
			items=items+[it]

	print('**************	List	*********************')
	out0 = open('Parse_List.txt','w')
	print('Items'+'\t'+'Items_URL'+'\t'+'Class'+'\t'+'Class_ID'+'\t'+'Items_ID'+'\t'+'Class_URL',file=out0)
	print([])
	out0.close()	
	out0 = open('Parse_List.txt','a')
	for i in range(0,len(items)):
		print(items[i]['category']+'\t'+items[i]['items_url']+'\t'+items[i]['_class']+'\t'+items[i]['class_id']+
			'\t'+items[i]['items_id']+'\t'+items[i]['class_url'],file=out0)
	out0.close()


