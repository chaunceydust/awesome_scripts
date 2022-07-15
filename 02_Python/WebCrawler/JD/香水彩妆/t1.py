#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import JDSpider
################	main	##########################

if __name__ == "__main__":
	key_word = ['book', 'e', 'channel', 'mvd', 'list']
	Base_url = 'https://list.jd.com'
	price_url = 'https://p.3.cn/prices/mgets?skuIds=J_'
	comment_url = 'https://club.jd.com/comment/productPageComments.action?productId=%s&score=0&sortType=5&page=%s&pageSize=10'
	favourable_url = 'https://cd.jd.com/promotion/v2?skuId=%s&area=1_72_2799_0&shopId=%s&venderId=%s&cat=%s'
	url='https://www.jd.com/allSort.aspx'
	key=['清洁用品','面部护肤','身体护理','口腔护理','女性护理','洗发护发','香水彩妆']
	logging.getLogger("requests").setLevel(logging.WARNING)  # 将requests的日志级别设成WARNING
############	抓取分类
	
	pc=JDSpider.JDSpider(url=url,key=key)
	cc=pc.parse_category()
	
	print('**************	Category	*********************')
	out0 = open('Parse_Category.txt','w')
	print('Items'+'\t'+'URL'+'\t'+'Class'+'\t'+'Class_ID',file=out0)
	out0.close()	
	out0 = open('Parse_Category.txt','a')
	for i in range(0,len(cc)):
		print(cc[i]['category']+'\t'+cc[i]['url']+'\t'+cc[i]['_class']+'\t'+cc[i]['class_id'],file=out0)
	out0.close()

