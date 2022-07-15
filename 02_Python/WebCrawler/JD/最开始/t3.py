#!/usr/bin/env python
# -*- coding: utf-8 -*-
import JDSpider
import Items
################	main	##########################

if __name__ == "__main__":
	key_word = ['book', 'e', 'channel', 'mvd', 'list']
	Base_url = 'https://list.jd.com'
	price_url = 'https://p.3.cn/prices/mgets?skuIds=J_'
	comment_url = 'https://club.jd.com/comment/productPageComments.action?productId=%s&score=0&sortType=5&page=%s&pageSize=10'
	favourable_url = 'https://cd.jd.com/promotion/v2?skuId=%s&area=1_72_2799_0&shopId=%s&venderId=%s&cat=%s'
	url='https://www.jd.com/allSort.aspx'
	key=['清洁用品','面部护肤','口腔护理','女性护理','洗发护发','香水彩妆']
###############
	items=[]
	files=open('Parse_List_new.txt')
	i=0
	for line in files.readlines():
		if i==0:
			i=2
			continue
		cutline=line.strip().split('\t')
		c0=Items.ProductsItem()
		c0['category']=cutline[0]
		c0['items_url']=cutline[1]
		c0['_class']=cutline[2]
		c0['class_id']=cutline[3]
		c0['items_id']=cutline[4]
#		c0['class_id']=cutline[5]
		items=items+[c0]
################ 抓取产品信息	
	out0 = open('Parse_Shop.txt','w')
	print('Shop_ID'+'\t'+'Shop_Name'+'\t'+'URL1'+'\t'+'URL2'+'\t'+'shopID'+'\t'+'venderID',file=out0)
	out0.close()	
	out0 = open('Parse_Shop.txt','a')
################################
	out1 = open('Parse_Product.txt','w')
	print('Name'+'\t'+'URL'+'\t'+'ProductID'+'\t'+'Category'+'\t'+'Items_id'+'\t'+'Class'+'\t'+'Class_id'+'\t'+'reallyPrice'+'\t'+'originalPrice'+
		'\t'+'description'+'\t'+'shopId'+'\t'+'favourableDesc1'+'\t'+'favourableDesc2'+'\t'+'goodRateShow'+'\t'+'poorRateShow'+'\t'+'generalRateShow'+
		'\t'+'goodCount'+'\t'+'generalCount'+'\t'+'poorCount'+'\t'+'commentCount'+'\t'+'defaultGoodCount'+'\t'+'averageScore',file=out1)
	out1.close()	
	out1 = open('Parse_Product.txt','a')
####################################
	out2 = open('Parse_Comment.txt','w')
	print('Name'+'\t'+'Count'+'\t'+'Category'+'\t'+'Class'+'\t'+'ProductID',file=out2)
	out2.close()	
	out2 = open('Parse_Comment.txt','a')
	for i in range(0,len(items)):
		n1=	items[i]['category']	
		u1=items[i]['items_url']
		c1=items[i]['_class']
		i1=items[i]['class_id']
		i2=items[i]['items_id']
		tt=JDSpider.JDSpider(url=u1,category=n1,class0=c1,class_id=i1,items_id=i2)
		print('Product:'+str(i)+' - '+str(len(items))+'\t\n'+u1)
		ffp=tt.parse_product()
		print(ffp[0]['_id']+'\t'+ffp[0]['name']+'\t'+ffp[0]['url1']+'\t'+ffp[0]['url2']+'\t'+ffp[0]['shopId']+'\t'+ffp[0]['venderId'],file=out0)
		print(ffp[1]['name']+'\t'+ffp[1]['items_url']+'\t'+ffp[1]['_id']+'\t'+ffp[1]['category']+'\t'+ffp[1]['items_id']+'\t'+ffp[1]['_class']+'\t'+ffp[1]['class_id']	
			+'\t'+ffp[1]['reallyPrice']+'\t'+ffp[1]['originalPrice']+'\t'+ffp[1]['description']+'\t'+ffp[1]['shopId']+'\t'+ffp[1]['favourableDesc1']
			+'\t'+ffp[1]['favourableDesc2']+'\t'+str(ffp[1]['goodRateShow'])+'\t'+str(ffp[1]['poorRateShow'])+'\t'+str(ffp[1]['generalRateShow'])
			+'\t'+str(ffp[1]['goodCount'])+'\t'+str(ffp[1]['generalCount'])+'\t'+str(ffp[1]['poorCount'])+'\t'+str(ffp[1]['commentCount'])
			+'\t'+str(ffp[1]['defaultGoodCount'])+'\t'+str(ffp[1]['averageScore']),file=out1)
		for cc in range(0,len(ffp[2])):
			print(ffp[2][cc]['name']+'\t'+str(ffp[2][cc]['count'])+'\t'+ffp[2][cc]['_class']+'\t'+
					ffp[2][cc]['category']+'\t'+str(ffp[2][cc]['productId']),file=out2)
	out0.close()
	out1.close()
	out2.close()


